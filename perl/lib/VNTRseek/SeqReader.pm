
=head1 NAME

VNTRseek::SeqReader - We use this module to process input read files
for TRF. The functions in this module take in file names and return
functions which act like reader objects, eg as in BioPerl.

=head1 SYNOPSIS

    use VNTRseek::SeqReader;

=head1 DESCRIPTION



=cut

package VNTRseek::SeqReader;
use v5.24;
use warnings;
use autodie;
use Try::Tiny;
use Carp;
use Exporter qw(import);
use POSIX qw(ceil);

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

    use File::Basename;

    # Get module's path safely (can't use FindBin reliably in a module)
    $self->{install_dir}
        = File::Basename::dirname( eval { ( caller() )[1] } ) . "/../../";
    no File::Basename;

    foreach my $param (
        qw( input_dir output_dir is_paired_end strip_454_TCAG
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
    # For new formats, simply add the name of the format here and require
    # that input file names have it as a prefix (eg, fasta files should
    # have prefix "fasta_")
    my @supported_format_names = qw(fasta fastq bam);
    my %supported_formats;
    @supported_formats{@supported_format_names} = @supported_format_names;

    # Dispatch table to call appropriate function for each format.
    # Uses the list @supported_format_names above for hash keys. When
    # adding a new format, add the name first to @supported_format_names,
    # and then add the reference to the reader function in the correct
    # place in the values list.
    my %reader_table;
    @reader_table{@supported_format_names}
        = ( \&read_fastaq, \&read_fastaq, \&read_bam );

    # TODO Simplify: we'll only support whatever seqtk supports and jusr
    # let seqtk handle decompression for us.
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
        if (@filenames = sort
            grep( /${supported_formats{$sf}}/ && -f "$_", @dircontents )
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
    if ( $input_format eq "bam" ) {
        @filenames = $self->init_bam( files => \@filenames );
        warn "BAM input. Will need to process "
            . scalar(@filenames)
            . " sets of reads from file.\n";
    }
    else {
        warn scalar(@filenames)
            . " supported files ($input_format format, $compression_msg) found in $self->{input_dir}\n";
    }

    # if ( $ENV{DEBUG} ) {
    #     use Data::Dumper;
    #     say Dumper( \@filenames );
    # }

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

    # if ( $ENV{DEBUG} ) {
    #     use Data::Dumper;
    #     say Dumper( \$next_input );
    # }
    if ($next_input) {
        return $self->{_reader}->(
            decom       => $self->{decom},
            install_dir => $self->{install_dir},
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
                if ( $self->{warn_454_TCAG} ) {
                    warn
                        "Read does not start with keyseq TCAG. Full sequence: $seq\n";
                }

                die
                    "Read does not start with keyseq TCAG. Full sequence: $seq\n";
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
    #     say "Queuing worker: $split_index, starting id: $start_id";
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

# DELETEME
sub wait_for_readers {
    my ( $self, %args ) = @_;

    $self->{read_loop}->await( $self->{read_future} );

    # $self->{redund_process}->close_when_empty();
    $self->{read_loop}->stop;
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
    my $self = shift;
    my %args = @_;

    my $read_href     = $args{input};
    my $output_prefix = $args{output_prefix};
    my $start_id      = $args{start_id};
    my $stat_diag     = "Writing with prefix = $output_prefix";

    use IPC::Run qw( run start pump finish timeout );

    my @trf_cmd        = ( $self->{trf_param}->@*, );
    my @trf2proclu_cmd = (
        $self->{trf2proclu_param}->@*,
        "-f", $start_id, "-o", "$output_prefix",
    );

    my @headers = ( keys $read_href->%* );

    my $trf_input = sub {
        state $i = 0;
        return if $i >= @headers;

        # warn "$headers[$i]\t" . $read_href->{$headers[$i]} . "\n";
        my $rchead
            = $headers[$i] . "_"
            . length( $read_href->{ $headers[$i] } )
            . "_RCYES";

        # NOTE: Must print reverse complement AFTER the forward read
        my $res
            = ">"
            . $headers[$i] . "\n"
            . $read_href->{ $headers[$i] }
            . "\n>$rchead\n"
            . reverse_complement( $read_href->{ $headers[$i] } );
        $i++;
        return $res;
    };

    my ( %reads, $trf_out );
    open my $trf_stdout, ">", "$output_prefix.leb36";
    my $proc_trf_output = sub {
        $trf_out .= $_[0];

        # Check if we got a whole line (or a chunk ending in a whole line)
        # That way, we can be sure we can get a whole header
        if ( $trf_out =~ /\n$/s ) {
            while ( $trf_out =~ m!^\d+\t([^/]+/[12])!mg ) {
                $reads{$1} = 1;

                # warn "Header " . $i++ . ": $1\n" if ($ENV{DEBUG});
            }
            print $trf_stdout $trf_out;
            $trf_out = "";
        }
    };

    my ( $trf_h, );

    try {
        $trf_h = start( \@trf_cmd, $trf_input, "|", \@trf2proclu_cmd,
            $proc_trf_output );
    }
    catch {
        # Try/catch just in case
        die "Cannot start TRF+trf2proclu pipe: $_\n";
    };

    pump $trf_h;
    finish $trf_h;

    my $reads_processed = keys $read_href->%*;

    if (%reads) {
        open my $reads_fh, ">", "$output_prefix.reads";

        for my $header ( keys %reads ) {
            unless ( exists $read_href->{$header} ) {
                die "Error: unlogged read (header $header) found in index. ",
                    "Possibly a bug in the sequence reading code. ",
                    "Please file a bug report.\n";
            }

            if ( $read_href->{$header} ) {
                say $reads_fh "$header\t" . $read_href->{$header};
                $read_href->{$header} = "";
            }
        }

        # .reads file gets one more line with the total number
        # of reads we read (different from reads with TRs count)
        say $reads_fh "totalreads\t$reads_processed";
        close $trf_stdout;
        close $reads_fh;
    }

    return { reads => $reads_processed, index => $args{index} };
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

    warn "Using $args{install_dir} for seqtk location.\n"
        if ( $ENV{DEBUG} );
    my $seqtk_bin = $args{install_dir} . "/seqtk";

    # Since we are using seqtk, use pipe open mode
    my $openmode = "-|";

    warn "Processing file " . $args{input} . "\n";
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

    # use Data::Dumper;
    # say "Input: " . Dumper($args{input});

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
