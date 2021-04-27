#!/usr/bin/env perl

# command line usage example:
#  ./produce_indist.pl full.leb36 filtered.leb36
# where 'full.leb36' is the full path to the unfiltered or larger
# reference set and 'filtered.leb36' is the full path to the filtered
# or smaller reference set.
#
# For determining indistinguishables, this script automates the process
# described in (TODO: cite), where multiple runs of profile clustering
# are performed with varying flank lengths (specifically 10, 20, and 50).
# Indistinguishables at each flank length are determined simply by
# checking if the filtered reference TR mapped to any other TR in the
# full, unfiltered set. The union of indistinguishables at all flank
# lengths is taken and written out to a file.
#
# This script allows configurable flank lengths (given as a comma
# separated list as an argument), and also collects data such as the
# number of TRs to which a TR was mapped. The latter is never used
# to make any determination, but may be necessary should the method be
# modified.
#

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use IO::Handle;
use FindBin;
use Cwd;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_ref_dbh make_refseq_db load_refprofiles_db run_redund);
use Data::Dumper;

my %opts = (
    perform_redund_full => 1,
    perform_redund_filt => 1,
    perform_psearch     => 1,
);
GetOptions( \%opts, "redo|r", );
die
    "Usage: $0 [options] full.leb36 filtered.leb36 flanklength1[,flanklength2,...] cutoff\n"
    . "\n"
    . "-r, --redo       Force recreating indist set and profile databases\n"
    unless @ARGV == 4;

my ( $full_set, $filtered_set, $flanklengths, $cutoff ) = @ARGV;
my $install_dir = $FindBin::Bin;
( my $filtered_basename = $filtered_set ) =~ s/.leb36$//;
( my $full_basename     = $full_set ) =~ s/.leb36$//;
my $redund_executable = "$install_dir/redund.exe";
my $proclu_executable = "$install_dir/psearch.exe";

#=<<Prepare SQLite dbs for reference set>>
#=<<First the filtered set>>
# Prepares the db with the ref set sequences, if not already done.
# Only need this for the filtered set.
warn "Initializing filtered set database (if needed).\n";
my $filtered_dbh = get_ref_dbh( $filtered_basename, { redo => $opts{redo} } );

#=<<Use databases to fetch rotindex files>>
# First the filtered set

# Get rotindex saved in db
my $get_rotindex = q{SELECT rotindex
    FROM files};

# For filtered set, fetch from the database, negate rid
my $get_profiles_sorted_q = q{SELECT -rid, length(pattern) AS patsize,
    printf("%.2f", copynum), proflen, proflenrc, profile, profilerc, nA, nC, nG, nT,
    printf("%s|%s", upper(substr(flankleft, -60)), upper(substr(flankright, 0, 61))) AS flanks
    FROM fasta_ref_reps JOIN ref_profiles USING (rid)
        JOIN minreporder USING (rid)
    ORDER BY minreporder.idx ASC};

# Create a temporary directory and write out filtered set
# and rotindex files there
my $tmpdir      = File::Temp->newdir();
my $tmpdir_name = $tmpdir->dirname;
my $tmp_filt_file
    = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );
my $tmp_full_file
    = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );

# Filtered file, ordered by min representation as given by redund
my $get_filtered_set_profiles_sth
    = $filtered_dbh->prepare($get_profiles_sorted_q);
$get_filtered_set_profiles_sth->execute;
open my $tmp_filt_file_fh, ">", $tmp_filt_file;
while ( my @fields = $get_filtered_set_profiles_sth->fetchrow_array ) {
    say $tmp_filt_file_fh join( " ", @fields );
}
close $tmp_filt_file_fh;

# Filtered rotindex
my ($rotindex_str) = $filtered_dbh->selectrow_array($get_rotindex);
die
    "Error getting rotindex file from filtered set. Try rerunning with --redo option\n"
    unless ($rotindex_str);
open my $tmp_rotindex_fh, ">", "${tmp_filt_file}.rotindex";

# Need to negate all indices
$rotindex_str =~ s/(\d+)/-$1/g;
print $tmp_rotindex_fh $rotindex_str;
close $tmp_rotindex_fh;

#=<<Then for the full file>>
my $tmp_full_dir;
warn "Initializing full set profiles (if needed).\n";
my $full_dbh = get_ref_dbh( $full_basename,
    { redo => $opts{redo}, skip_refseq => 1 } );

# Full file, ordered by min representation as given by redund
open my $full_set_fh, "<", $full_set;
my %full_hash;
while ( my $line = <$full_set_fh> ) {
    my ( $rid, @fields ) = split /\s+/, $line;
    $full_hash{$rid} = \@fields;
}
close $full_set_fh;

$get_profiles_sorted_q = q{SELECT rid
    FROM ref_profiles
        JOIN minreporder USING (rid)
    ORDER BY minreporder.idx ASC};
my $get_full_set_profiles_sth = $full_dbh->prepare($get_profiles_sorted_q);
$get_full_set_profiles_sth->execute;
open my $tmp_full_file_fh, ">", $tmp_full_file;
while ( my ($rid) = $get_full_set_profiles_sth->fetchrow_array ) {
    say $tmp_full_file_fh join( " ", $rid, @{ $full_hash{$rid} } );
}
close $tmp_full_file_fh;
%full_hash = ();

# Full rotindex
($rotindex_str) = $full_dbh->selectrow_array($get_rotindex);
die
    "Error getting rotindex file from filtered set. Try rerunning with --redo option\n"
    unless ($rotindex_str);
my $tmp_full_rotindex = "${tmp_full_file}.rotindex";
open $tmp_rotindex_fh, ">", $tmp_full_rotindex;

# No need to negate all indices
print $tmp_rotindex_fh $rotindex_str;
close $tmp_rotindex_fh;

# print "Press enter to continue...";
# my $dummy = <STDIN>;

#=<<Now, run psearch.exe for each flanklength given>>
my @flens = split /,/, $flanklengths;
my ( %indist, %files_seen_in );
for my $fl (@flens) {
    my $out_prefix = "reads${fl}_F_${cutoff}";
    my $proclu_cmd
        = qq($proclu_executable $tmp_filt_file $tmp_full_file $install_dir/eucledian.dst $cutoff 0 0 -r $fl -m);
    warn $proclu_cmd . "\n";
    system($proclu_cmd);
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    # Process each map file, get all the TRs in the filtered set which mapped
    # to a TR other than itself. Filtered set TRs are negative.
    open my $mapping, "<", "$tmp_full_file.map";
    my %num_links_per_tr;
    while ( my $line = <$mapping> ) {
        $line =~ /(-\d+)(['"])=>(\d+):.*/;

        # if ( ( exists( $trlist{$1} ) ) && ( -$1 != $3 ) ) {
        if ( -$1 != $3 ) {
            if ( $ENV{DEBUG} ) {
                warn
                    "$1 mapped to TR other than itself at flank length $fl: $3\n";
            }

            # Check if seen as indist for this param set. Only
            # increment global counter if we haven't seen a mapping
            # to another TR for this parameter set yet. Keeps track
            # of how many flank lengths agreed this was indist.
            # (Counts are not currently used.)
            $indist{$1}++ unless ( exists $num_links_per_tr{$1} );

            # Mark that we've seen as indist this param set and count
            # how many map links seen so far per TR.
            $num_links_per_tr{$1}++;
        }
    }
}

# Find all elements of %indist where the count matches the length of the
# @flens array. These are the indists.
$filtered_dbh->begin_work();
$filtered_dbh->do(
    q{CREATE TEMPORARY TABLE indist (
    rid integer PRIMARY KEY)}
) or die "Couldn't do statement: $DBI::errstr\n";
my $insert_indist_sth = $filtered_dbh->prepare(
    q{INSERT INTO indist
        (rid) VALUES (?)}
);
$filtered_dbh->commit;

# NB Can get rid of this after we eliminate the need for indist files.
my $bname = basename( $filtered_set, ".leb36" );
open my $indist_out, ">", getcwd . "/$bname.indist";
for my $tr ( sort keys(%indist) ) {
    say $indist_out $tr;
    $insert_indist_sth->execute( -$tr );
}
close $indist_out;

# End NB

warn "Updating indistinguishable TRs...\n";
$filtered_dbh->begin_work;
my $reset_indist = q{UPDATE fasta_ref_reps
    SET is_singleton = 1, is_dist = 1, is_indist = 0};
my $reset_count = $filtered_dbh->do($reset_indist);
unless ($reset_count) {
    warn "0 or undef TRs reset when updating indistinguishables. "
        . "This could mean there was a problem, so check your database "
        . "when this script exits using the query: "
        . "'SELECT is_indist, COUNT(*) FROM fasta_ref_reps GROUP BY is_indist;' "
        . "from the sqlite3 command line to check that you have TRs in both "
        . "categories. Contact us if you need help.\n";
}
else {
    warn "$reset_count TRs reset to default indist status.\n";
}

my $update_indist = q{UPDATE fasta_ref_reps
    SET is_singleton = 0, is_dist = 0, is_indist = 1,
    WHERE EXISTS (
    SELECT rid FROM indist t2
    WHERE fasta_ref_reps.rid = t2.rid
)};
my $indist_count = $filtered_dbh->do($update_indist);
warn "$indist_count indistinguishable TRs...\n";
$filtered_dbh->commit;
$filtered_dbh->disconnect;

