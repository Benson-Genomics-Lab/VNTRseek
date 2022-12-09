#!/usr/bin/env perl

# command line usage example:
#  ./run_trf_ng.pl fasta outputdir
# where 6 is the number of files to process in one batch,
#   fasta is the input directory containing zipped files,
#   and outputdir is the directory for output

use v5.24;
use strict;
use warnings;
use IO::Handle;
use Parallel::ForkManager;
use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config set_statistics);
use SeqReader;

# Arguments
die "Usage: run_trf_ng.pl expects 8 arguments.\n"
    unless scalar @ARGV >= 8;

my $input_dir        = $ARGV[0];
my $output_dir       = $ARGV[1];
my $cnf              = $ARGV[2];
my $trf_param        = $ARGV[3];
my $trf2proclu_param = $ARGV[4];
my $strip_454_TCAG   = $ARGV[5];
my $is_paired_end    = $ARGV[6];
my $max_processes    = $ARGV[7];


# Database info for set_statistics and SeqReader init
my %run_conf = get_config("CONFIG", $cnf);


$trf_param =~ s/['"]([^'"]+)['"] //;
my $trf_bin = $1;
$trf2proclu_param =~ s/['"]([^'"]+)['"] //;
my $trf2proclu_bin = $1;


my $seq_reader = SeqReader->new(
    seqtk            => "$FindBin::RealBin/seqtk",
    input_dir        => $input_dir,
    output_dir       => $output_dir,
    is_paired_end    => $is_paired_end,
    reads_split      => 1e6,
    trf_param        => [ $trf_bin, split /\s+/, $trf_param ],
    trf2proclu_param => [ $trf2proclu_bin, split /\s+/, $trf2proclu_param ],
);


my %trf_res = (
    reads             => 0,
    num_trs_ge7       => 0,
    num_trs           => 0,
    num_reads_trs_ge7 => 0,
    num_reads_trs     => 0,
);

my %timetrf = ();


print "Will use $max_processes processes.\n";

my $pm = Parallel::ForkManager->new($max_processes);
$pm->run_on_finish(
    sub {
        my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $res ) = @_;
        if ($res) {
            for my $k ( keys $res->%* ) {
                if (exists $trf_res{$k}) { $trf_res{$k} += $res->{$k}; }
            }

            ( $ENV{DEBUG} )
                && warn "Process $ident read $res->{reads} reads.\n";
        }
    }
);
$pm->set_waitpid_blocking_sleep(0);

# Asynchronously collect input from input reader
# As soon as we have read_split reads, launch a new TRF process,
# or wait if not enough workers available.
my $split_index = 0;
READS:
while ( my $reads = $seq_reader->get_reads() ) {
    $pm->start($split_index) and ( $split_index++, next READS );

    # Child code
    $pm->finish( 0, { reads => 0 } ) unless $reads->%*;

    my $start_id = ( $split_index * $seq_reader->{reads_split} ) + 1;
    my $res = $seq_reader->run_trf(
        output_prefix => "$seq_reader->{output_dir}/$split_index",
        index         => $split_index,
        input         => $reads,
        start_id      => $start_id,
    );
    $pm->finish( 0, $res );
}

print "Finished reading. Waiting for TRF processes.\n";
$pm->wait_all_children();


print "Processing complete -- processed $split_index file(s)"
    . " and $trf_res{reads} reads.\n";

set_statistics({
    NUMBER_READS             => $trf_res{'reads'},
    NUMBER_TRS_IN_READS_GE7  => $trf_res{'num_trs_ge7'},
    NUMBER_TRS_IN_READS      => $trf_res{'num_trs'},
    NUMBER_READS_WITHTRS_GE7 => $trf_res{'num_reads_trs_ge7'},
    NUMBER_READS_WITHTRS     => $trf_res{'num_reads_trs'},
});
