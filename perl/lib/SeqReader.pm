=head1 NAME

SeqReader - We use this module to process input read files
for TRF. The functions in this module take in file names and return
functions which act like reader objects, eg as in BioPerl.

=head1 SYNOPSIS

    use SeqReader;

=head1 DESCRIPTION

=cut

package SeqReader;
use v5.24;
use warnings;
use autodie;
use Try::Tiny;
use Carp;
use Exporter "import";
use POSIX "ceil";

# use IO::Async::Function;
# use IO::Async::Stream;
# use IO::Async::Process;
# use IO::Async::Loop;

=item C<new>

Initializes a new instance of the sequence reader

## @brief      { function_description }
##
## @param      my    { parameter_description }
##
## @return     { description_of_the_return_value }
##

=cut

##
sub new {
    my ( $class, %attrs ) = @_;
    my $self = bless {}, $class;

    # In case I want to do something with $self first
    $self->_init(%attrs);

    return $self;
}

sub _init {
    my ( $self, %params ) = @_;

    foreach my $param (
        qw( seqtk input_dir output_dir is_paired_end strip_454_TCAG
        warn_454_TCAG reads_split trf_param trf2proclu_param )
        )
    {
        $self->{$param} = delete $params{$param};
    }

    $self->{reads_split} = 1e6 unless defined $self->{reads_split};

    $self->_init_input_list(%params);
}

sub _init_input_list {
    my ( $self, %params ) = @_;

    # get a list of input files
    croak "Input directory does not exist" unless -d $self->{input_dir};
    opendir( my $dirhandle, $self->{input_dir} );
    my @dircontents = readdir($dirhandle);
    closedir($dirhandle);
    croak "No input files found" unless @dircontents;

    # Prepend the directory to the returned file names
    @dircontents = map { $self->{input_dir} . "/$_" } @dircontents;

    ( $ENV{DEBUG} )
        && warn join( ",", @dircontents ) . "\n";
    my ( @filenames, $input_format );
    my $compression = '';

    # List all supported file extensions and input formats here. Order of
    # input formats is in priority order: first format if found is used.
    # For new formats, simply add the name of the format here and its extensions.
    my @supported_format_names = qw(fasta fastq sam bam cram);
    my %supported_formats;
    @supported_formats{@supported_format_names} = (
        qr/fa|fasta/, qr/fq|fastq/,
        qr/sam/, qr/bam/, qr/cram/
    );


    # Dispatch table to call appropriate function for each format.
    # Uses the list @supported_format_names above for hash keys. When
    # adding a new format, add the name first to @supported_format_names,
    # and then add the reference (\&) to the reader function in the correct
    # place in the values list.
    my %reader_table;
    @reader_table{@supported_format_names}
        = ( \&read_fastaq, \&read_fastaq, \&read_bam, \&read_bam, \&read_bam);

    # Similarly for compression formats, add a new one to the list, then
    # add the regex and command needed to extract in the correct position
    # in the values list
    my @compressed_format_names = qw(targz gzip tarbzip bzip tarxz xz);
    my %compressed_formats;
    @compressed_formats{@compressed_format_names} = (
        qr/tgz|tar\.gz/,               qr/gz/,
        qr/tbz|tbz2|tar\.bz|tar\.bz2/, qr/bz|bz2/,
        qr/txz|tar\.xz/,               qr/xz/
    );

    # All decompression commands end with a space, for convenience
    my %decompress_cmds;
    @decompress_cmds{@compressed_format_names} = (
        "tar xzfmO ",
        "gunzip -c ",
        "tar xjfmO ",
        "bzip2 -c ",
        "tar xJfmO ",
        "xzcat "
    );
    while ( my $sf = shift @supported_format_names ) {
        if ( @filenames = sort
             grep ( /(?:${supported_formats{$sf}})(?:\.|$)/ && -f "$_", @dircontents )
           )
        {
            $input_format = $sf;

            # Remove bam.bai files
            @filenames = grep ( !/\.bai$/, @filenames );

            # Determine compression format
            while ( my ( $cf, $cf_re ) = each %compressed_formats ) {
                if ( $filenames[0] =~ /.*\.(?:${cf_re})/ ) {
                    $compression = $cf;
                    last;
                }
            }
            last;
        }
    }

    my $compression_msg
        = ($compression)
        ? "compressed as $compression"
        : "assuming uncompressed";
    croak "No supported files found in $self->{input_dir}. Exiting\n"
        if ( @filenames == 0 );

    # If BAM, init files list
    if ($input_format eq "sam" or
        $input_format eq "bam" or
        $input_format eq "cram" ) {

        @filenames = $self->init_bam( files => \@filenames );
        print "SAM/BAM/CRAM input. Will need to process "
            . scalar(@filenames)
            . " sets of reads from file.\n";
    }
    else {
        print scalar(@filenames)
            . " supported files ($input_format format, $compression_msg) found in $self->{input_dir}\n";
    }

    $self->{num_inputs}   = scalar @filenames;
    $self->{inputs}       = \@filenames;
    $self->{input_format} = $input_format;
    $self->{decom}   = ($compression) ? $decompress_cmds{$compression} : '';
    $self->{_reader} = $reader_table{$input_format};

    # Init first reader
    $self->{reader} = $self->_next_reader;

    # $self->{input_index} = 0;
}

=item I<_next_reader()>

...

=cut

sub _next_reader {
    my $self       = shift;
    my $file_index = $self->{num_inputs} - @{ $self->{inputs} };
    my $next_input = shift( @{ $self->{inputs} } );

    if ($next_input) {
        return $self->{_reader}->(
            decom       => $self->{decom},
            seqtk       => $self->{seqtk},
            input       => $next_input,
            file_index  => $file_index,
        );
    }

    return;
}

=item I<get_reads()>

...

=cut

sub get_reads {
    my $self          = shift;
    my $rc_read_count = 0;

    ( $ENV{DEBUG} ) && warn "Num inputs $self->{num_inputs}\n";
    my %read_hash;
    my $read_count;

    # my $future;
    # my $split_index = 0;

    # As long as there are readers to get input from
    while ( defined $self->{reader}
        || ( $self->{reader} = $self->_next_reader ) )
    {

        # While reader returns sequence records, process and accumulate.
        while ( my ( $header, $seq ) = $self->{reader}->() ) {
            if ( exists $read_hash{$header} ) {
                die "Error: duplicate reads in input. Make sure your"
                    . "input only consists of unique reads. This can happen"
                    . "if your input contains alternative alignments of the"
                    . "same sequence after converting from BAM/CRAM.\n";
            }

            # Trim tags, if needed, and prepare reverse complement
            if ( $self->{strip_454_TCAG} && ( $seq !~ s/^TCAG//i ) ) {
                die "Read $header does not start with keyseq TCAG. Full sequence: $seq\n";
            }

            $read_hash{$header} = $seq;

            # Once we've accumulated reads_split records, queue
            # an instance of TRF function with the records read and
            # clear records list.
            if ( ( ++$read_count % $self->{reads_split} ) == 0 ) {
                ( $ENV{DEBUG} )
                    && warn "Return reads\n";
                return \%read_hash;
            }
        }

        # If exited loop, finished with this reader
        if ( $ENV{DEBUG} ) {
            warn "Finished a reader loop\n";
        }
        $self->{reader} = undef;
    }

    # Return any reads (if read less than read_split reads)
    if (%read_hash) {
        ( $ENV{DEBUG} )
            && warn "Return reads\n";
        return \%read_hash;
    }

    return;

    # if (@read_list) {
    #     my $start_id = ( $split_index * $self->{reads_split} ) + 1;
    #     #print "Queuing worker: $split_index, starting id: $start_id\n";
    #     my $f = $read_function->call(
    #         args => [
    #             {   output_prefix    => "$self->{output_dir}/$split_index",
    #                 index            => $split_index,
    #                 trf_param        => $self->{trf_param},
    #                 trf2proclu_param => $self->{trf2proclu_param},
    #                 input            => [@read_list],
    #                 start_id         => $start_id,
    #             }
    #         ],
    #         )->on_done(
    #         sub {
    #             my $res = shift;
    #             my $split_index = $res->{index};
    #             ( $ENV{DEBUG} )
    #                 && warn
    #                 "Process $split_index read $res->{reads} reads.\n";
    #         }
    #         );
    #     $split_index++;
    #     @read_list = ();
    #     $future = $f unless $future;
    # }

    # ( $ENV{DEBUG} )
    #     && warn "Reading process finished. ",
    #     "Waiting for TRF processes.", "\n";
    # $self->{read_future} = $future;

    # $self->{read_loop}->await_all(@futures);
}

=item I<run_trf()>

Takes input/output dirs, the TRF and TRF2proclu command calls with
parameters, flags for reversing reads and 454 tag handling, file
format and compression format names, the number of files processed
so far, and a list of names of files in that format. Then
returns a function which, when called, returns a stream to sequences
which can be used to run TRF/TRF2PROCLU.

Relies on a global, reader_table, which links sequence format names
to functions which can read those formats. Functions in that table
must return a sub which returns a list of exactly two values each
time it is called, and undef when there are no more sequences. The
two values are a FASTA header and a sequence string, in that order.

=cut

sub run_trf {
    use IPC::Run qw( run new_chunker );

    my $self           = shift;
    my %args           = @_;
    my $read_href      = $args{input};
    my $output_prefix  = $args{output_prefix};
    my $start_id       = $args{start_id};
    my @trf_cmd        = ( $self->{trf_param}->@*, );
    my @trf2proclu_cmd = (
        $self->{trf2proclu_param}->@*,
        "-f", $start_id, "-o", "$output_prefix",
    );

    my @headers   = ( keys $read_href->%* );
    my $num_reads = 0;

    my $trf_input = sub {
        return if $num_reads >= @headers;

        my $rchead
            = $headers[$num_reads] . "_"
            . length( $read_href->{ $headers[$num_reads] } )
            . "_RCYES";

        # NOTE: Must print reverse complement AFTER the forward read
        my $res
            = ">"
            . $headers[$num_reads] . "\n"
            . $read_href->{ $headers[$num_reads] }
            . "\n>$rchead\n"
            . ( reverse($read_href->{ $headers[$num_reads] })
                    =~ tr/ACGTacgt/TGCAtgca/r );
        $num_reads++;
        return $res;
    };

    my ( %tr_reads, $trf_h );
    open my $trf_stdout, ">", "$output_prefix.index";

    my $proc_trf_output = sub {

        chomp( my $line = $_[0] );
        print $trf_stdout $_[0];
        warn "Chunk: $line\n" if ( $ENV{DEBUG} );

        my @f = split /\t/, $line;
        if ( @f == 7 ) {
            $tr_reads{ $f[1] } = 1;

            warn "Header: $f[1]\n" if ( $ENV{DEBUG} );
            return;
        }

        croak "Invalid output from TRF2PROCLU ($_[0])\n";
    };

    try {
        run(\@trf_cmd, $trf_input, "|", \@trf2proclu_cmd,
            '>', new_chunker, $proc_trf_output );
    }
    catch {
        # Try/catch just in case
        die "Cannot start TRF+trf2proclu pipe: $_\n";
    };

    close $trf_stdout;

    my %res = ( reads => $num_reads, index => $args{index} );

    # If there were any reads with TRs, dump those reads to a file
    if (%tr_reads) {
        open my $reads_fh, ">", "$output_prefix.reads";
        for my $header ( keys %tr_reads ) {
            if ( $read_href->{$header} ) {
                print $reads_fh "$header\t" . $read_href->{$header} . "\n";
                $read_href->{$header} = "";
            }
        }
        close $reads_fh;
    }


    # Get stats from indexhist file
    open my $fh, "<", "$output_prefix.indexhist";
    chomp( my $line = <$fh> );
    @res{qw(
        num_trs_ge7
        num_trs
        num_reads_trs_ge7
        num_reads_trs
        )} = split(/\t/, $line);
    close $fh;

    return \%res;
}

=item I<reverse_complement()>

Takes a DNA sequence string and returns the reverse complement.

=cut

sub reverse_complement {
    return reverse( $_[0] ) =~ tr/ACGTacgt/TGCAtgca/r;
}

=item I<read_fastaq()>

FASTA/Q file reader. Requires seqtk.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This particular sub does NOT modify the $files_to_process value.

=cut

sub read_fastaq {
    my %args = @_;

    warn "Using $args{seqtk} for seqtk location.\n"
        if ( $ENV{DEBUG} );
    my $seqtk_bin = $args{seqtk};

    # Since we are using seqtk, use pipe open mode
    my $openmode = "-|";

    # this will be a problem given many input fastas
    print "Processing file " . $args{input} . "\n";
    my $filename = '"' . $args{input} . '"';

    if ( $args{decom} ) {
        $filename = $args{decom} . $filename . "| $seqtk_bin seq -a -S";
    }
    else {
        $filename = "$seqtk_bin seq -a -S " . $filename;
    }

    # warn "Filename/command = '$filename'\n";
    open my ($fasta_fh), $openmode, $filename;

   # Consume first empty record because of the way $/ splits the FASTA format.
   # <$fasta_fh>;
    return sub {
        local $/ = "\n>";
        my $fasta_rec = <$fasta_fh>;
        return () unless ($fasta_rec);
        my ( $header, $seq ) = split( /\n+/, $fasta_rec );
        chomp $header;
        $header =~ s/^>//;
        $header =~ s/\s+$//;
        chomp $seq;

        # Add a number (the file index) to the read
        # if the read information cannot be determined
        # from the read header
        state $need_idx = !(

            # Illumina BaseSpace FASTQ header
            ( $header =~ / [12]:[YN]:\d+:(\d+|[ACGTacgt]+)/ )
            ||

            # Other flag seen to indicate pair
            ( $header =~ /\/[12]/ )
        );
        ($need_idx) && ( $header .= " vs=$args{file_index}" );

        # warn "header: '$header'";
        # my $seq = join( "", @seqlines );

        # warn "seq: '$seq'";
        return ( $header, $seq );
    };
}

=item I<init_bam()>

Initialize arrayref filelist for BAM file reader.

Given the input dir, a flag indicating if the input is paired-end
reads, and a list of all files, return a list of generated samtools
commands for each portion on the input BAM file(s).

The returned list can be used to replace filelist in the caller,
which will then be passed to the rest of the program for running TRF.

Requires samtools as an external dependency.

=cut

sub init_bam {
    my ( $self, %args ) = @_;
    my @samcmds;

    # $start is always 1
    # TODO Use bedtools bamtofastq?
    # For samtools, can use "-L <bedfile> -M" and have
    # separate command for all unmapped
    my $samviewcmd   = "samtools view ";
    my $unpairedflag = "-F 1"
        ; # Probably not required: only single-end fragments in a single-end template anyway
    my $firstsegflag      = "-f 64";
    my $lastsegflag       = "-f 128";
    my $unmappedflag      = "-f 4";
    my $badmapflag        = "-F 256 -F 2048";
    my $unmapped_template = "*";

    for my $file ( @{ $args{files} } ) {
        my $bamfile = $file;

        my $samviewflags = $badmapflag;
        if ( $self->{is_paired_end} ) {
            my $cmd = join( ' ',
                $samviewcmd, $samviewflags, $firstsegflag, $bamfile );
            push @samcmds, { cmd => $cmd, pair => "/1" };
            $cmd = join( ' ',
                $samviewcmd, $samviewflags, $lastsegflag, $bamfile );
            push @samcmds, { cmd => $cmd, pair => "/2" };
        }
        else {
            my $cmd = join( ' ',
                $samviewcmd, $samviewflags, $unpairedflag, $bamfile );
            push @samcmds, { cmd => $cmd, pair => "" };
        }
    }

    return @samcmds;
}

=item I<read_bam()>

BAM file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

Requires samtools as an external dependency.

=cut

sub read_bam {
    my %args = @_;

    # warn "$current_idx\n";
    my $cmdhash = $args{input};
    ( $ENV{DEBUG} )
        && warn "Processing bam chunk using: " . $cmdhash->{cmd} . "\n";

    local $SIG{PIPE} = sub { die "Error in samtools pipe: $?\n" };
    open my $samout, "-|", $cmdhash->{cmd};
    return sub {
        my ( $header, $lpos, $seq ) = ( "", -1, "" );

        # # Skip redundant reads from previous regions. Start with
        # # lpos = -1 so that we always read at least one line.
        # # We also must enter the loop if processing the unmapped
        # # reads.
        # while ( ( $cmdhash->{start} eq "unmapped" )
        #     || $lpos < $cmdhash->{start} )
        # {
        my $bam_rec = <$samout>;
        return () unless ($bam_rec);

        (   $header, undef, undef, $lpos, undef,
            undef,   undef, undef, undef, $seq
        ) = split( /\t/, $bam_rec );

        #     # If we got a read, and this is an unmapped read, just
        #     # jump out of loop.
        #     last if ( $cmdhash->{start} eq "unmapped" );
        # }
        return ( "$header" . $cmdhash->{pair}, $seq );
    };
}

=item I<read_FORMAT() stub>

Example file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This is an example of how to write a sub which reads a given file
in a certain format and returns a subroutine which works as expected
by fork_proc().

sub read_FORMAT {
    my ( $input_dir, $compression, $current_file, $files_to_process, $filelist )
        = @_;

    # Here make some decisions. For instance, if we expect to process multiple
    # files, such as FASTA or FASTQ files and unlike BAM files, then choose
    # the file to read using $current_file. Otherwise, logically split
    # the input and decide which split this call will work on. You may have to
    # write the next section before this one in the case of multiple files.

    # In the case of files like BAM files, where there is only one input
    # file but the file is indexed in regions which can be read directly,
    # save the number of splits needed into the variable pointed to by
    # $files_to_process. NOT needed if there are multiple input files,
    # one for each child process.

    # Open a file handle to the file or pipe we will read from.
    # If file might be compressed, use the $compression value and
    # the decompress_cmds hash to construct the input pipe command.
    my $input_fh = open...

    # Handle any errors, and die if there is a fatal problem so this file
    # or split isn't considered further and another process can jump in on
    # the next item in the list.

    # Return an anonymous subroutine which will use that file handle.
    return sub {
        # Read from input file handle
        my $record = <$input_fh>;
        # Return empty list if nothing else to read
        return () unless ($record);

        # Process $record as needed, read more lines, etc (see other functions
        # for examples).

        # Return a list consisting of the header, prefixed by a ">", and the
        # sequence string.
        return(">" . $header, $seq);
    }
}

=back

=head1 LICENSE

This is released under the Artistic 
License. See L<perlartistic>.

=head1 AUTHOR

Yozen Hernandez
Yevgeniy Gelfand

=head1 SEE ALSO

L<perl>

=cut

1;
