#!/usr/bin/env perl

# command line usage example:
#  ./run_trf.pl 6 fasta
# where 6 is the number of files to process in one batch
# and fasta is the input directory containing zipped files
#

use strict;
use warnings;
use feature 'state';
use Getopt::Std;
use IO::Handle;
use Carp;
use FindBin;

# this is where the pipeline is installed
use lib "$FindBin::RealBin/lib";    # must be same as install dir!
use ProcInputReads
    qw(fork_proc init_bam formats_regexs compressed_formats_regexs set_install_dir);
our $install_dir = "$FindBin::RealBin";
set_install_dir($install_dir);

my $files_processed  = 0;    # files processed
my $files_to_process = 0;    # Either: the actual number of files to process
     # OR the number of splits of a BAM file (and others?)
my %p;    # associates forked pids with output pipe pids

my %opts;
getopts( 'p:t:u:swr', \%opts );

if ( !scalar(@ARGV) || !defined( $opts{'t'} ) || !defined( $opts{'u'} ) ) {
    die "Usage:
        perl $0 [-p <processes>] -t <trf> -u <trf2proclu> -s (strip_454_TCAG) -w (warn_454_TCAG) -r (IS_PAIRED_READS) <input_dir> <output_dir>\n";
}

my $TRF_PARAM        = $opts{'t'};
my $TRF2PROCLU_PARAM = $opts{'u'};
my $strip_454_TCAG   = ( exists $opts{'s'} && $opts{'s'} ) ? 1 : 0;
my $warn_454_TCAG    = ( exists $opts{'w'} && $opts{'w'} ) ? 1 : 0;
my $IS_PAIRED_READS  = ( exists $opts{'r'} && $opts{'r'} ) ? 1 : 0;
my $max_processes    = ( exists $opts{'p'} && $opts{'p'} ) ? $opts{'p'} : 2;

my $reverse_read = 1; # if 1, each read will be reversed and processed as well

my $HEADER_SUFFIX = ""
    ;   # set by program automatically for paired reads to distinguish headers


# Input directory
my $input_dir = $ARGV[0];
die "Need to provide input directory\n" if !defined $input_dir;
die "Input directory does not exist\n"  if !-d $input_dir;

# Output directory
my $output_dir = $ARGV[1];
die "Need to provide output directory\n" if !defined $output_dir;

#die "Output directory already exists -- please delete it before running this program\n" if -e $output_dir;
mkdir $output_dir;

# get a list of input files
opendir( my $dirhandle, $input_dir );
my @dircontents = readdir($dirhandle);
closedir($dirhandle);

if ( $ENV{DEBUG} ) {
    warn join( ",", @dircontents ) . "\n";
}

# Get all supported files. See note above on priority of input formats
my @filenames;
my ( $input_format, $compression );

# Determine sequence format
my %supported_formats_regexs = formats_regexs();
while ( my ( $sf, $pat_re ) = each %supported_formats_regexs ) {
    if (@filenames = sort
        grep( /${pat_re}/ && -f "$input_dir/$_", @dircontents )
        )
    {
        $input_format = $sf;

        # Remove bam.bai files
        @filenames = grep ( !/\.bai$/, @filenames );

        # Determine compression format
        my %cmp_formats_regexs = compressed_formats_regexs();
        while ( my ( $cf, $cf_re ) = each %cmp_formats_regexs ) {
            if ( $filenames[0] =~ /.*\.(?:${cf_re})/ ) {
                $compression = $cf;
                last;
            }
        }
        last;
    }
}

$files_to_process = @filenames;
my $compression_msg
    = ($compression)
    ? "compressed as $compression"
    : "assuming uncompressed";
die "0 supported files found in $input_dir. Exiting\n"
    if $files_to_process == 0;
warn
    "$files_to_process supported files ($input_format format, $compression_msg) found in $input_dir\n";

# /proc/cpuinfo is found on Linux/cygwin only: use other methods for other platforms
if ( $max_processes == 0 ) {
    open my $cpuinfo, "<", "/proc/cpuinfo";
    if ($cpuinfo) {
        $max_processes = scalar( map /^processor/, <$cpuinfo> );
        warn "$max_processes CPU core(s) detected\n";
        if ( $max_processes == 0 ) {
            warn "Unknown formatting in /proc/cpuinfo: Assuming single CPU\n";
            $max_processes = 1;
        }
    }
    else {
        warn "Could not open /proc/cpuinfo: Assuming single CPU\n";
        $max_processes = 1;
    }
}

# If BAM, init files list
if ( $input_format eq "bam" ) {
    @filenames = init_bam( $input_dir, $IS_PAIRED_READS, \@filenames );
    $files_to_process = @filenames;
    warn
        "BAM input. Will need to process $files_to_process sets of reads from file.\n";
}

$max_processes
    = ( $files_to_process < $max_processes )
    ? $files_to_process
    : $max_processes;
warn "Will use $max_processes processes\n";

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $max_processes; $i++ ) {
    $p{ fork_proc(
            $input_dir,        $output_dir,        $TRF_PARAM,
            $TRF2PROCLU_PARAM, $reverse_read,      $strip_454_TCAG,
            $warn_454_TCAG,    $input_format,      $compression,
            $files_processed,  \$files_to_process, \@filenames
        )
    } = 1;
    $files_processed++;
}

my $num_files = $files_to_process;

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

        # Only spawn more processes if there are still more files/splits
        # to process.
        if ( $files_processed < $num_files ) {
            $p{ fork_proc(
                    $input_dir,         $output_dir,
                    $TRF_PARAM,         $TRF2PROCLU_PARAM,
                    $reverse_read,      $strip_454_TCAG,
                    $warn_454_TCAG,     $input_format,
                    $compression,       $files_processed,
                    \$files_to_process, \@filenames
                )
            } = 1;
            $files_processed++;
        }
    }
    else {
        die
            "ERROR: Process wait()ed on was not found in our process list. PID=$pid\n";
    }
}

warn "Processing complete -- processed $files_processed file(s).\n";
my $exstring = qq(find "${output_dir}" -type f -size 0 -delete);
system($exstring);
