
=head1 NAME

ProcInputReads - We use this module to process input read files for
TRF. The functions in this module take in file names and return
functions which act like reader objects, eg as in BioPerl.

=head1 SYNOPSIS

    use ProcInputReads;
    my 

=head1 DESCRIPTION

This module does not really exist, it
was made for the sole purpose of
demonstrating how POD works.

=cut

package ProcInputReads;
use strict;
use warnings;
use Carp;
use feature 'say';
use Exporter qw(import);
use POSIX qw(ceil);

our @EXPORT_OK
    = qw(fork_proc init_bam get_reader formats_regexs compressed_formats_regexs set_install_dir);

my $install_dir;

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

my $records_before_split = 1e6;

=head2 Methods

=over 12

=item C<set_install_dir>

Returns a hash of sequence format names and the regular expressions
used to determine it from a file name.

=cut

sub set_install_dir {
    $install_dir = shift;
    croak
        "Argument to set_install_dir must be a defined and non-empty path to the install directory."
        unless ( ( defined $install_dir ) && ( $install_dir ne "" ) );
    warn "Set install_dir to $install_dir"
        if ( $ENV{DEBUG} );
}

=item C<formats_regexs>

Returns a hash of sequence format names and the regular expressions
used to determine it from a file name.

=cut

sub formats_regexs {
    return %supported_formats;
}

=item C<check_trf_pipe_close>

Checks the return status of the trf+trf2proclu pipe, and handles
each case appropriately.

=cut

#** @function private check_trf_pipe_close ($errno, $child_err, $output_prefix, $current_file, $current_fragment)
# @brief Checks the return status of the trf+trf2proclu pipe, and handles
# each case appropriately.
#
# Wraps the perl error checking convention for calls to system
# for our trf+trf2proclu pipe, and handles the possible exit states
# of the command.
# @params errno required integer - $!, the perl variable set after a failure by a system call
# @params child_err required integer - $?, the perl variable for an error in a child process
# @params output_prefix required string - A string indicating base name of a file produced by the pipe
# @params current_file requred integer - Index of the currently running process
# @params current_fragment required integer - Index of the current fragment of the input for this process
# @retval integer or FALSE - Exits on error, return value of "unlink" on success.
#*
sub check_trf_pipe_close {
    my ( $errno, $child_err, $output_prefix, $current_file,
        $current_fragment )
        = @_;
    my $trf_pipe_es  = $child_err >> 8;
    my $trf_pipe_sig = $child_err & 127;

    # Here the process trf2proclu finished with non-0 AND there was some other
    # problem, since $! is not 0.
    if ($!) {
        warn "Error closing trf+trf2proclu pipe: $!\n";
        exit(255);
    }

    # TODO Fix: error 14 might be returned because of some error in TRF
    elsif ( $trf_pipe_es == 1 || $trf_pipe_es == 14 ) {

        # No TRs processed, remove files
        warn
            "Did not find any TRs (process $current_file, part $current_fragment). Removing leb36 file...\n";
        unlink("$output_prefix.leb36");
    }

    # Here trf2proclu exited with non-0 status due to some error
    elsif ( $trf_pipe_es > 1 ) {
        warn "trf+trf2proclu pipe has returned $trf_pipe_es\n";
        unlink("$output_prefix.leb36");
        exit($trf_pipe_es);
    }
}

=item C<write_reads>

Write out reads with TRs to files in the output directory.

=cut

##** @function write_reads ($output_prefix, $reads_processed, $read_hash)
##
## @brief      Writes reads with TRs to output directory. File names
##             have same prefix as the corresponding index file.
##
## @param      $output_prefix    The output file prefix
## @param      $reads_processed  The number of reads processed before
##             calling this function by this process.
## @param      $read_hash        A hashref with read headers as keys,
##             and sequences as values.
##
## @return     none. Dies on error.
##
sub write_reads {
    my ( $output_prefix, $reads_processed, $read_hash ) = @_;
    open my $index_fh, "<", "$output_prefix.index";
    open my $reads_fh, ">", "$output_prefix.reads";
    while ( my $line = <$index_fh> ) {
        my @fields = split /\t/, $line;
        my $head = $fields[1];
        if ( exists $read_hash->{$head}
            && ( $read_hash->{$head}->{written} == 0 ) )
        {
            print $reads_fh $head . "\t" . $read_hash->{$head}->{seq} . "\n";
            $read_hash->{$head}->{written} = 1;
        }
        elsif ( !exists $read_hash->{$head} ) {
            die "Error: unlogged read found in index."
                . "Possibly a bug in the sequence reading"
                . "code. Please file a bug report.\n";
        }
    }

    # .reads file gets one more line with the total number
    # of reads we read (different from reads with TRs count)
    say $reads_fh "totalreads\t$reads_processed";
    close $index_fh;
    close $reads_fh;
}

=item C<compressed_formats_regexs>

Returns a hash of sequence format names and the regular expressions
used to determine it from a file name.

=cut

sub compressed_formats_regexs {
    return %compressed_formats;
}

=item I<fork_proc()>

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

sub fork_proc {
    my ($input_dir,        $output_dir,       $trf_param,
        $trf2proclu_param, $reverse_read,     $strip_454_TCAG,
        $warn_454_TCAG,    $format,           $compression,
        $current_file,     $files_to_process, $filelist
    ) = @_;

    defined( my $pid = fork() )
        or die "Unable to fork: $!\n";
    if ( $pid == 0 ) {    #Child

        my $reader = $reader_table{$format}->(
            $input_dir, $compression, $current_file,
            $files_to_process, $filelist
        );
        exit unless $reader;

        my $current_fragment = 0;
        my $output_prefix = "$output_dir/${current_file}_$current_fragment";
        my $stat_diag     = "Running child, current_file = $current_file";
        if ( ref( $filelist->[$current_file] ) eq "HASH" ) {
            while ( my ( $k, $v ) = each %{ $filelist->[$current_file] } ) {
                $stat_diag .= qq{, $k = "$v"};
            }
        }
        else {
            $stat_diag .= " ($filelist->[$current_file])";
        }

        warn $stat_diag, "...\n";

        # TODO Error checking if TRF, in the start of the pipe, breaks down
        local $SIG{PIPE} = sub { die "Error in trf+trf2proclu pipe: $?\n" };
        open my $trf_pipe,
            "|$trf_param | $trf2proclu_param -o '$output_prefix.index' > '$output_prefix.leb36'"
            or die "Cannot start TRF+trf2proclu pipe: $!\n";

        # TODO Need way of logging TRF output?
        # open my $logfile, ">", "$output_prefix.log"
        #     or die "Error opening file $output_prefix.log: $!\n";
        # $logfile->autoflush;
        # my $debug_reads_processed = 0;

        my $reads_processed = 0;
        my %reads;
        while ( my ( $header, $body ) = $reader->() ) {
            if ( $reads_processed && ( $reads_processed++ % $records_before_split ) == 0 ) {
                if ( !close $trf_pipe ) {
                    check_trf_pipe_close( $!, $?, $output_prefix,
                        $current_file, $current_fragment );
                }

                # scan output file and write out reads with TRs to new file
                if ( -e "$output_prefix.index" && -s "$output_prefix.index" )
                {
                    write_reads( $output_prefix, $reads_processed, \%reads );
                }

                # Reset reads hash and reads_processed. Prevents
                # memory bloat in the former case, and lets us record
                # only the reads processed this split in the latter
                %reads           = ();
                $reads_processed = 0;
                $output_prefix
                    = "$output_dir/${current_file}_" . ++$current_fragment;
                open $trf_pipe,
                    "|$trf_param | $trf2proclu_param -o "
                    . "'$output_prefix.index' > '$output_prefix.leb36'"
                    or die "Cannot start TRF+trf2proclu pipe: $!\n";
            }
            if ( exists $reads{$header} ) {
                die "Error: duplicate reads in input. Make sure your"
                    . "input only consists of unique reads. This can happen"
                    . "if your input contains alternative alignments of the"
                    . "same sequence after converting from BAM/CRAM.\n";
            }
            else {
                $reads{$header} = { seq => $body, written => 0 };
            }

            # warn "header: $header\nbody: $body\n"
            #     if ($ENV{DEBUG});

            # say $logfile $data[0] . "\n" . $data[1];

            pipe_to_trf( $reverse_read, $strip_454_TCAG, $warn_454_TCAG,
                $trf_pipe, ">" . $header, $body );

            # $trf_pipe, $header, $body, $debug_reads_processed );
            # $debug_reads_processed++;
        }

        # Normally, close() returns false for failure of a pipe.
        if ( !close $trf_pipe ) {
            check_trf_pipe_close( $!, $?, $output_prefix, $current_file,
                $current_fragment );
        }

        # scan output file and write out reads with TRs to new file
        if ( -e "$output_prefix.index" && -s "$output_prefix.index" ) {
            write_reads( $output_prefix, $reads_processed, \%reads );
        }

        exit;
    }
    else {    # Parent
              # warn "Running parent, out_counter = $out_counter...\n";
        return $pid;
    }
}

=item I<get_reader()>

Takes input dir, file format and compression format names,
the number of files processed so far, and a list of names of files
in that format. Then returns a function which, when called, returns
a stream to sequences from the read file.

=cut

sub get_reader {
    my ( $input_dir, $format, $compression, $current_file,
        $files_to_process, $filelist )
        = @_;

    croak
        "Calling script must call set_install_dir with the install dir path before using other functions in this package\n"
        unless ( ( defined $install_dir ) && ( $install_dir ne "" ) );

    my $reader = $reader_table{$format}->(
        $input_dir, $compression, $current_file, $files_to_process, $filelist
    );

    exit unless $reader;
    return $reader;
}

=item I<reverse_complement()>

Takes a DNA sequence string and returns the reverse complement.

=cut

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

=item I<pipe_to_trf()>

Takes flags related to reversing reads and handling 454 tags, a file
handle to a TRF process pipe to TRF2proclu, a FASTA header, and a
FASTA sequence (plain sequence string) as input.

Processes input FASTA records to reverse them or remove 454 tags,
as dictated by user options. Prints processed output directly to
given file handle to TRF pipeline.

=cut

sub pipe_to_trf {
    my ( $reverse_read, $strip_454_TCAG, $warn_454_TCAG, $trf_fh, $header,
        $body )

        # $body, $reads_processed )
        = @_;

    # warn "Processing header $header";

    if ( $strip_454_TCAG && ( $body !~ s/^TCAG//i ) ) {
        if ($warn_454_TCAG) {
            warn
                "Read does not start with keyseq TCAG. Full sequence: $body\n";
        }
        else {
            die
                "Read does not start with keyseq TCAG. Full sequence: $body\n";
        }
    }

    # NOTE: Must print reverse complement AFTER the forward read
    say $trf_fh "$header\n$body";

    if ( $reverse_read && $header ne "" ) {

        # if ( ($reads_processed > 1) && $reverse_read && $header ne "" ) {
        say $trf_fh $header . "_"
            . length($body)
            . "_RCYES\n"
            . reverse_complement($body);
    }
}

=item I<read_fasta()>

FASTA file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This particular sub does NOT modify the $files_to_process value.

=cut

sub read_fasta {
    my ( $input_dir, $compression, $current_file,
        $files_to_process, $filelist )
        = @_;

# warn "Need to process " . scalar(@$filelist) . " files, working on $current_file\n";

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    warn "Processing file " . $filelist->[$current_file] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$current_file] . '"';

    # $current_file contains how many files processed so far.
    # Use to index into filelist
    # warn $filename;
    local $/ = ">";

    # warn "Filename/command = '$filename'\n";
    open my $fasta_fh, $openmode, $filename
        or die "Error opening file " . $filename;

   # Consume first empty record because of the way $/ splits the FASTA format.
    <$fasta_fh>;
    return sub {
        local $/ = ">";
        my $fasta_rec = <$fasta_fh>;
        return () unless ($fasta_rec);
        my ( $header, @seqlines ) = split( /\n+/, $fasta_rec );
        chomp $header;
        chomp @seqlines;

        # warn "header: '$header'";
        my $seq = join( "", @seqlines );

        # warn "seq: '$seq'";
        return ( ">" . $header, $seq );
    };
}

=item I<read_fastq()>

FASTQ file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This particular sub does NOT modify the $files_to_process value.

=cut

sub read_fastq {

    # Code modified from https://www.biostars.org/p/11599/#11657
    my ( $input_dir, $compression, $current_file,
        $files_to_process, $filelist )
        = @_;

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    warn "Processing file " . $filelist->[$current_file] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$current_file] . '"';

    # $current_file contains how many files processed so far.
    # Use to index into filelist
    # warn "Filename/command = '$filename'\n";
    open my $fastq_fh, $openmode, $filename
        or die "Error opening file " . $filename;

    my $aux         = undef;
    my $seq_counter = 0;

    return sub {

        # warn "Seq number: " . $seq_counter++ . "\n";
        $aux = [ undef, 0 ] unless ( defined $aux );
        return () if ( $aux->[1] );
        if ( !defined( $aux->[0] ) ) {
            while (<$fastq_fh>) {
                chomp;
                if ( substr( $_, 0, 1 ) eq '>' || substr( $_, 0, 1 ) eq '@' )
                {
                    $aux->[0] = $_;
                    last;
                }
            }
            if ( !defined( $aux->[0] ) ) {
                $aux->[1] = 1;
                return ();
            }
        }

        # warn "After header: $.\n";
        my $name = $aux->[0] =~ /^.(\S+)/ ? $1 : '';
        my $seq = '';
        my $c;
        $aux->[0] = undef;
        while (<$fastq_fh>) {
            chomp;
            $c = substr( $_, 0, 1 );
            last if ( $c eq '>' || $c eq '@' || $c eq '+' );
            $seq .= $_;
        }

        # warn "After seq: $.\n";
        $aux->[0] = $_;
        $aux->[1] = 1 if ( !defined( $aux->[0] ) );
        croak "Next line is not qual header: c = $c" if ( $c ne '+' );

        # return (">" . $name, $seq) if ($c ne '+');
        my $qual = '';
        while (<$fastq_fh>) {
            chomp;
            $qual .= $_;
            if ( length($qual) >= length($seq) ) {
                $aux->[0] = undef;
                last;
            }
        }

        # warn "After qual: $.\n";
        $aux->[1] = 0;
        return ( ">" . $name, $seq );
    };
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
    my ( $input_dir, $compression, $current_file,
        $files_to_process, $filelist )
        = @_;

    unless ( $current_file < $$files_to_process ) {
        return undef;
    }

    warn "Using $install_dir for seqtk location.\n"
        if ( $ENV{DEBUG} );
    my $seqtk_bin = $install_dir . "/seqtk";

# warn "Need to process " . scalar(@$filelist) . " files, working on $current_file\n";

    # Since we are using seqtk, use pipe open mode
    my $openmode = "-|";

    warn "Processing file " . $filelist->[$current_file] . "\n";
    my $filename = '"' . "$input_dir/" . $filelist->[$current_file] . '"';

    if ($compression) {
        $filename
            = $decompress_cmds{$compression}
            . $filename
            . "| $seqtk_bin seq -a -S";
    }
    else {
        $filename = "$seqtk_bin seq -a -S " . $filename;
    }

    # $current_file contains how many files processed so far.
    # Use to index into filelist
    # warn $filename;
    local $/ = ">";

    # warn "Filename/command = '$filename'\n";
    open my $fasta_fh, $openmode, $filename
        or die "Error opening file " . $filename;

   # Consume first empty record because of the way $/ splits the FASTA format.
    <$fasta_fh>;
    use feature 'state';
    return sub {
        local $/ = ">";
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
        ($need_idx) && ( $header .= " vs=$current_file" );

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
    my ( $input_dir, $is_paired_end, $filelist ) = @_;
    my @samcmds;

    # $start is always 1
    my $samviewcmd   = "samtools view ";
    my $unpairedflag = "-F 1"
        ; # Probably not required: only single-end fragments in a single-end template anyway
    my $firstsegflag      = "-f 64";
    my $lastsegflag       = "-f 128";
    my $unmappedflag      = "-f 4";
    my $badmapflag        = "-F 256 -F 2048";
    my $unmapped_template = "*";

    for my $file (@$filelist) {
        next if ( $file =~ /\.bai$/ );
        my $bamfile = "$input_dir/" . $file;

        # warn "$bamfile\n";
        # my $baifile = $bamfile . ".bai";

        # warn "$baifile\n";

        # Requires samtools to be installed/available
        # Get all regions in the bam file
        # Don't save sequence with 0 reads
        my @regions = grep ( !/0\s+0$/, qx( samtools idxstats "$bamfile" ) );
        if ( $? == -1 ) {
            croak "failed to execute: $!\n";
        }
        elsif ( $? & 127 ) {
            croak sprintf(
                "command died while processing file '%s' with signal %d, %s coredump. Does a bai file exist?\n",
                $bamfile,
                ( $? & 127 ),
                ( $? & 128 ) ? 'with' : 'without'
            );
        }

        # Don't do more if we've exhausted the file list
        # warn "Regions: " . scalar(@regions) . "\n";

# We need to read the regions in the bam file and construct samtools commands
# These are all saved in an array which are interated through like the FASTA files before.
# Then we catch the output of these through a file handle, and process into FASTA.
# TODO Change how we segment these: should we divide the sequences into fixed read portions? Divide each chromosome into even parts? Etc..
        for my $r (@regions) {
            my ( $chr, $end, $num_aln, $num_unaln ) = split /\s+/, $r;
            my $unmapped = ( $chr eq $unmapped_template );

            # next if ( ( $num_aln + $num_unaln ) == 0 );

            my $samviewflags = ($unmapped) ? $unmappedflag : $badmapflag;

            # Redundant code, simplfiy later
            if ($unmapped) {
                if ($is_paired_end) {
                    my $cmd = join( ' ',
                        $samviewcmd, $samviewflags, $firstsegflag, $bamfile,
                        "'*'" );
                    push @samcmds,
                        { cmd => $cmd, pair => "/1", start => "unmapped" };
                    $cmd = join( ' ',
                        $samviewcmd, $samviewflags, $lastsegflag, $bamfile,
                        "'*'" );
                    push @samcmds,
                        { cmd => $cmd, pair => "/2", start => "unmapped" };
                }
                else {
                    my $cmd = join( ' ',
                        $samviewcmd, $samviewflags, $unpairedflag, $bamfile,
                        "'*'" );
                    push @samcmds,
                        { cmd => $cmd, pair => "", start => "unmapped" };
                }
            }
            else {
                my @end_coords
                    = map { $_ * 10e6 } ( 1 .. ceil( $end / 10e6 ) );

                for my $i ( 0 .. $#end_coords ) {
                    my $start_coord = ( 1 + ( $i * 10e6 ) );
                    my $region
                        = "$chr:" . $start_coord . "-" . $end_coords[$i];

                    if ($is_paired_end) {
                        my $cmd = join( ' ',
                            $samviewcmd, $samviewflags, $firstsegflag,
                            $bamfile, $region );
                        push @samcmds,
                            {
                            cmd   => $cmd,
                            pair  => "/1",
                            start => $start_coord
                            };
                        $cmd = join( ' ',
                            $samviewcmd, $samviewflags, $lastsegflag,
                            $bamfile, $region );
                        push @samcmds,
                            {
                            cmd   => $cmd,
                            pair  => "/2",
                            start => $start_coord
                            };
                    }
                    else {
                        my $cmd = join( ' ',
                            $samviewcmd, $samviewflags, $unpairedflag,
                            $bamfile, $region );
                        push @samcmds,
                            {
                            cmd   => $cmd,
                            pair  => "",
                            start => $start_coord
                            };
                    }
                }
            }
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
    my ( $input_dir, $compression, $current_idx, $files_to_process, $cmdlist )
        = @_;

# Return undef if the index "$current_idx" exceeds the number of samtools commands
    return undef
        unless ( $current_idx < @$cmdlist );

    # warn "$current_idx\n";
    my $cmdhash = $cmdlist->[$current_idx];
    ( $ENV{DEBUG} )
        && warn "Processing bam chunk using: " . $cmdhash->{cmd} . "\n";

    local $SIG{PIPE} = sub { die "Error in samtools pipe: $?\n" };
    open my $samout, "-|", $cmdhash->{cmd}
        or die "Error opening samtools pipe: $!\n";
    return sub {
        my ( $header, $lpos, $seq ) = ( "", -1, "" );

        # Skip redundant reads from previous regions. Start with
        # lpos = -1 so that we always read at least one line.
        # We also must enter the loop if processing the unmapped
        # reads.
        while ( ( $cmdhash->{start} eq "unmapped" )
            || $lpos < $cmdhash->{start} )
        {
            my $bam_rec = <$samout>;
            return () unless ($bam_rec);

            (   $header, undef, undef, $lpos, undef,
                undef,   undef, undef, undef, $seq
            ) = split( /\t/, $bam_rec );

            # If we got a read, and this is an unmapped read, just
            # jump out of loop.
            last if ( $cmdhash->{start} eq "unmapped" );
        }
        return ( $header . $cmdhash->{pair}, $seq );
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
