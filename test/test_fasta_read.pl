#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Carp;

# Test reading FASTA files and outputting a sequence and header.
die "Usage: $0 <format> <compression> <fasta file 1> ... <fasta file N>"
    unless @ARGV >= 3;
my $max_processes = 2;
my $reverse_read = 1; # if 1, each read will be reversed and processed as well
my $strip_454_TCAG = 0;
my $warn_454_TCAG  = 0;
my $TRF_EXE        = $ENV{HOME} . "/src/VNTRseek/build/trf409-ngs.linux.exe";
my $TRF_PARAM      = "'$TRF_EXE' - 2 5 7 80 10 50 2000 -d -h -ngs";
my $TRF2PROCLU_EXE
    = $ENV{HOME} . '/src/VNTRseek/build/src/trf2proclu-ngs/trf2proclu-ngs.exe';
my $TRF2PROCLU_PARAM
    = "'$TRF2PROCLU_EXE' -f 1 -m 2 -s 5 -i 7 -p 7 -l 5";
my ( $format, $compression, @filenames ) = @ARGV;

# For this testing script only
$compression = ( $compression eq "none" ) ? undef : $compression;
my %decompress_cmds = (
    targz   => "tar xzfmO ",
    gzip    => "gunzip -c ",
    tarbzip => "tar xjfmO ",
    bzip    => "bzip2 -c ",
    tarxz   => "tar xJfmO ",
    xz      => "xzcat "
);
my %reader_table = (
    fasta => \&read_fasta,
    fastq => \&read_fastq,
    bam   => \&read_bam
);

# Setting here because needed in full program
my $input_dir   = "";
my $output_dir  = ".";
my $out_counter = 0;
my %p;

for ( my $i = 0; $i < $max_processes; $i++ ) {
    $p{ fork_proc( $format, $compression, $out_counter, \@filenames ) } = 1;
    $out_counter++;
}

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {

    # check return value
    my ( $rc, $sig, $core ) = ( $? >> 8, $? & 127, $? & 128 );
    if ($core) {
        warn "run_trf process $pid dumped core\n";
        exit(1000);
    }
    elsif ( $sig == 9 ) {
        warn "run_trf process $pid was murdered!\n";
        exit(1001);
    }
    elsif ( $rc != 0 ) {
        warn "run_trf process $pid has returned $rc!\n";
        exit($rc);
    }

    if ( $p{$pid} ) {

        # one instance has finished processing -- start a new one
        delete $p{$pid};

     # For BAM files (and maybe other formats?) the parent does not know
     # how many total processes are needed to run the whole file. For these
     # formats, reader functions need to write out a temporary file signalling
     # the last processes that needs to run has already begun.
        if ( -e "$output_dir/trf_alldone" ) {
            unlink("$output_dir/trf_alldone");
            last;
        }
        else {
            $p{ fork_proc( $format, $compression, $out_counter, \@filenames ) } = 1;
            $out_counter++;
        }
    }
    else {
        die "ERROR: Do not remember process PID=$pid\n";
    }
}

# while ( my $fasta_reader = read_bam( "", \@filenames ) ) {
#     last if $out_counter == 4;

#     # my $fasta_reader = read_fasta( 8, "", \@filenames );
#     open my $outfile, ">", "$output_dir/$out_counter.out"
#         or die "Error opening file $output_dir/$out_counter.out: $!\n";
#     while ( my @data = $fasta_reader->() ) {
#         say $outfile $data[0] . "\n" . $data[1];
#     }
#     say "";
#     $out_counter++;
# }

sub fork_proc {
    my ( $format, $compression, $files_processed, $filelist ) = @_;
    defined( my $pid = fork() )
        or die "Unable to fork: $!\n";
    if ( $pid == 0 ) {    #Child

        my $reader
            = $reader_table{$format}
            ->( $compression, $files_processed, $filelist );
        my $output_prefix = "$output_dir/$files_processed";
        exit unless $reader;
        warn "Running child, files_processed = $files_processed...\n";
        # TODO Error checking if TRF, in the start of the pipe, breaks down
        local $SIG{PIPE} = sub { die "Error in trf+trf2proclu pipe: $?\n" };
        open my $trf_pipe,
            "|$TRF_PARAM | $TRF2PROCLU_PARAM -o '$output_prefix.index' > '$output_prefix.leb36'"
            or die "Cannot start TRF+trf2proclu pipe: $!\n";
        # TODO Need way of logging TRF output?
        # open my $logfile, ">", "$output_prefix.log"
        #     or die "Error opening file $output_prefix.log: $!\n";
        # $logfile->autoflush;
        my $debug_trs_found = 0;

        while ( my ($header, $body) = $reader->() ) {
            # say $logfile $data[0] . "\n" . $data[1];
            pipe_to_trf( $trf_pipe, $header, $body );
        }

        # Normally, close() returns false for failure of a pipe. If the only problem
        # was that the exit status of the pipe was non-0, then $! == 0.
        # Important because trf2proclu returns non-0 on success.
        if ( !close $trf_pipe ) {
            # Here the process trf2proclu finished with non-0 AND there was some other
            # problem, since $! is not 0.
            if ($!) {
                warn "Error closing trf+trf2proclu pipe: $!\n";
                exit(1002);
            }
            # Here trf2proclu exited with non-0 status but that was the only issue,
            # so just report that value.
            elsif ( $? < -2 ) {
                warn "trf+trf2proclu pipe has returned $?\n";
                exit($?);
            }
            # TODO Doesn't account for TRF's negative return values
        }

        # Check exit error


        exit;
    }
    else {    # Parent
              # warn "Running parent, out_counter = $out_counter...\n";
        return $pid;
    }
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

# Processes input FASTA records to reverse them or remove 454 tags
# Prints processed output directly to given file handle to TRF
# process.
sub pipe_to_trf {
    my ( $trf_fh, $header, $body ) = @_;
    # warn "Processing header $header";

    if ( $reverse_read && $header ne "" ) {
        say $trf_fh $header . "_" . length($body) . "_RCYES";
        say $trf_fh reverse_complement($body);
    }

    # FASTA header
    say $trf_fh $header;
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
    say $trf_fh $body;
}

sub read_fasta {
    my ( $compression, $files_processed, $filelist ) = @_;
    # warn "Need to process " . scalar(@$filelist) . " files, working on $files_processed\n";

    # Don't do more if we've exhausted the file list
    if ( $files_processed >= @$filelist ) {
        system("touch '$output_dir/trf_alldone'");
        return undef;
    }
    
    # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    warn "Processing file " . $filelist->[$files_processed] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$files_processed] . '"';

    # $files_processed contains how many files processed so far.
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

sub read_fastq {
    # Code modified from https://www.biostars.org/p/11599/#11657
    my ( $compression, $files_processed, $filelist ) = @_;

    # Don't do more if we've exhausted the file list
    if ( $files_processed >= @$filelist ) {
        system("touch '$output_dir/trf_alldone'");
        return undef;
    }
    
    # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    warn "Processing file " . $filelist->[$files_processed] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$files_processed] . '"';

    # $files_processed contains how many files processed so far.
    # Use to index into filelist
    # warn "Filename/command = '$filename'\n";
    open my $fastq_fh, $openmode, $filename
        or die "Error opening file " . $filename;

    my $aux = undef;
    my $seq_counter = 0;

    return sub {
        # warn "Seq number: " . $seq_counter++ . "\n";
        $aux = [undef, 0] unless (defined $aux);
        return () if ($aux->[1]);
        if (!defined($aux->[0])) {
            while (<$fastq_fh>) {
                chomp;
                if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                    $aux->[0] = $_;
                    last;
                }
            }
            if (!defined($aux->[0])) {
                $aux->[1] = 1;
                return ();
            }
        }
        # warn "After header: $.\n";
        my $name = $aux->[0] =~ /^.(\S+)/? $1 : '';
        my $seq = '';
        my $c;
        $aux->[0] = undef;
        while (<$fastq_fh>) {
            chomp;
            $c = substr($_, 0, 1);
            last if ($c eq '>' || $c eq '@' || $c eq '+');
            $seq .= $_;
        }
        # warn "After seq: $.\n";
        $aux->[0] = $_;
        $aux->[1] = 1 if (!defined($aux->[0]));
        croak "Next line is not qual header: c = $c" if ($c ne '+');
        # return (">" . $name, $seq) if ($c ne '+');
        my $qual = '';
        while (<$fastq_fh>) {
            chomp;
            $qual .= $_;
            if (length($qual) >= length($seq)) {
                $aux->[0] = undef;
                last;
            }
        }
        # warn "After qual: $.\n";
        $aux->[1] = 0;
        return ( ">" . $name, $seq );
    };
}

# Requires samtools
sub read_bam {
    my ( $compression, $files_processed, $filelist ) = @_;
    my $bamfile = "$input_dir/" . $filelist->[0];

    # warn "$bamfile\n";
    my $baifile = $bamfile . ".bai";

    # warn "$baifile\n";
    my $unmapped_template = "*";
    my @samcmds;

    # Check if .bai file exists and then run MakeBedFiles.jar
    die
        "Error reading bam file: corresponding bai file required for processing bam files."
        unless ( -e -r $baifile );

    # Requires samtools to be installed/available
    # Get all regions in the bam file
    my @regions = qx(
            samtools idxstats "$bamfile"
        );

    # Don't do more if we've exhausted the file list
    # warn "Regions: " . scalar(@regions) . "\n";

# We need to read the regions in the bam file and construct samtools commands
# These are all saved in an array which are interated through like the FASTA files before.
# Then we catch the output of these through a file handle, and process into FASTA.
    for my $r (@regions) {
        my ( $chr, $end, $num_aln, $num_unaln ) = split /\s+/, $r;
        my $unmapped = ( $chr eq $unmapped_template );
        # Don't save sequence with 0 reads
        next if ( ( $num_aln + $num_unaln ) == 0 );

        # $start is always 1
        my $region = "$chr:1-$end";
        my $scmd
            = "samtools view "
            . ( ($unmapped) ? "-f 4 " : "" )
            . $bamfile
            . ( ($unmapped) ? "" : " $region" );
        push @samcmds, $scmd;
    }

    if ( $files_processed >= @samcmds ) {
        system("touch '$output_dir/trf_alldone'");
        return undef;
    }
    # warn "$files_processed\n";
    warn "Processing bam chunk using: " . $samcmds[$files_processed] . "\n";

    local $SIG{PIPE} = sub { die "Error in samtools pipe: $?\n" };
    open my $samout, "-|", $samcmds[ $files_processed++ ]
        or die "Error opening samtools pipe: $!\n";
    return sub {
        my $bam_rec = <$samout>;
        return () unless ($bam_rec);

        # print ">" $1 "\n" $10
        my ($header, undef, undef, undef, undef,
            undef,   undef, undef, undef, $seq
        ) = split( /\s+/, $bam_rec );
        return ( ">" . $header, $seq );
        }
}
