#!/usr/bin/env perl

# sets a number of db stat variables
#
# command line usage example:
#  ./setdbstats.pl reads_profiles_folder reference_folder reads_profile_folder_clean dbname dblogin dbpass dbhost
# where inputfile is the main cluster file
#

use strict;
use warnings;
use Cwd;

use FindBin;
use File::Basename;

use List::Util qw[min max];

use lib "$FindBin::RealBin/lib";    # must be same as install dir!

use vutil qw( get_config get_ref_dbh set_statistics get_statistics );

my $argc = @ARGV;
if ( $argc < 4 ) {
    die
        "Usage: setdbstats.pl reads_profiles_folder reads_profile_folder_clean dbsuffix run_dir\n";
}

my $readpf   = $ARGV[0];
my $rpfc     = $ARGV[1];
my $DBSUFFIX = $ARGV[2];
my $run_dir  = $ARGV[3];

####################################

my %stats;
my $exstring;
my $input;
my $rc;

my %run_conf = get_config( $DBSUFFIX, $run_dir );
my $dbh = get_ref_dbh($run_conf{REFERENCE}, { readonly => 1 } );

( $stats{NUMBER_REF_TRS} )
    = $dbh->selectrow_array(q{SELECT COUNT(*) FROM fasta_ref_reps});

# open($input, "-|", "wc -l $reffolder/reference.leb36.rotindex | tail -1");
# $rc = <$input>;
# if ($rc =~ /(\d+)/) {
#   $stats{NUMBER_REFS_TRS_AFTER_REDUND} = $1;
# }
# close($input);

( $stats{NUMBER_REFS_TRS_AFTER_REDUND} )
    = $dbh->selectrow_array(
    q{SELECT COUNT(*) FROM ref_profiles WHERE redund = 0});

$dbh->disconnect;

open( $input, "-|", "cat $rpfc/*.rotindex | wc -l" );
$rc = <$input>;
if ( $rc =~ /(\d+)/ ) {
    $stats{NUMBER_TRS_IN_READS_AFTER_REDUND}      = $1;
    $stats{NUMBER_READS_WITHTRS_GE7_AFTER_REDUND} = $1;
}
close($input);

# Get these stats. If none are set, or at least one is unset,
# re-read the index files and set.
# Its recommended that copies of databases are made by simply
# copying the run directory and database, so somehow skipping
# the step where these are set is not expected to be an issue.
my $num_tr_stats = get_statistics(
    qw(NUMBER_TRS_IN_READS
        NUMBER_TRS_IN_READS_GE7
        NUMBER_READS_WITHTRS_GE7
        NUMBER_READS_WITHTRS)
);
if ( !ref $num_tr_stats || grep { !defined $_ || $_ == 0 }
    values $num_tr_stats->%* )
{
    $stats{NUMBER_TRS_IN_READS}      = 0;
    $stats{NUMBER_TRS_IN_READS_GE7}  = 0;
    $stats{NUMBER_READS_WITHTRS_GE7} = 0;
    $stats{NUMBER_READS_WITHTRS}     = 0;

    opendir( my $dirhandle, "$readpf" );
    my @dircontents = grep( /\.indexhist$/, readdir($dirhandle) );
    closedir($dirhandle);
    for my $f (@dircontents) {
        open my $fh, "<", "$readpf/$f";
        my $line = <$fh>;
        chomp $line;
        my @fields = split /\t/, $line;
        $stats{NUMBER_TRS_IN_READS_GE7}  += $fields[0];
        $stats{NUMBER_TRS_IN_READS}      += $fields[1];
        $stats{NUMBER_READS_WITHTRS_GE7} += $fields[2];
        $stats{NUMBER_READS_WITHTRS}     += $fields[3];
        close $fh;
    }
}

# open( $input, "-|", "cat $readpf/*.indexhist | wc -l" );
# $rc = <$input>;
# if ( $rc =~ /(\d+)/ ) {
#     $stats{NUMBER_TRS_IN_READS} = $1;
# }
# close($input);

# $rc = qx(./ge7.pl $readpf/*.index 2>/dev/null);
# if ( $rc =~ /(\d+) (\d+) (\d+)/ ) {
#     $stats{NUMBER_TRS_IN_READS_GE7}  = $1;
#     $stats{NUMBER_READS_WITHTRS_GE7} = $2;
#     $stats{NUMBER_READS_WITHTRS}     = $3;
# }
# close($input);

# open( $input, "-|", "./ge7.pl $readpf/*.indexhist" );
# $rc = <$input>;
# if ( $rc =~ /(\d+) (\d+) (\d+)/ ) {
#     $stats{NUMBER_READS_WITHTRS} = $3;
# }
# close($input);

# open( $input, "-|", "cat $rpfc/*.rotindex | wc -l" );
# $rc = <$input>;
# if ( $rc =~ /^(\d+)/ ) {
#     $stats{NUMBER_READS_WITHTRS_GE7_AFTER_REDUND} = $1;
# }
# close($input);

# Get config for run and save stats
set_statistics( \%stats );

1;
