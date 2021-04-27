#!/usr/bin/env perl

# command line usage example:
#  ./remove_sw_fpvntrs.pl <TODO>
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
use feature qw(say state);
use Getopt::Std;
use IO::Handle;
use FindBin;
use Cwd;
use DBI;

# use File::Copy;
use File::Basename;

# use File::Temp qw/ tempfile tempdir /;
use lib "$FindBin::RealBin/lib";
use vutil qw( get_config get_dbh );

# Get input files
die "Usage: $0 <dbsuffix> <db backend> [new reference set name]"
    unless ( @ARGV && ( @ARGV <= 3 ) );
my ( $dbsuffix, $backend, $refset_suffix ) = @ARGV;
my $run_dir = getcwd();

# Connect to DB
my %run_conf = get_config( $dbsuffix, $run_dir );
my ( $login, $pass, $host ) = @run_conf{qw(LOGIN PASS HOST)};
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

# Retrieve from the database the list of all reads mapping to a reftr
# my $all_reads_mapped_query = q{SELECT DISTINCT map.refid AS reftrid, SUBSTRING_INDEX(fasta_reads.head, "_", -1) AS origintrid
#   FROM map INNER JOIN replnk ON replnk.rid=map.readid
#       INNER JOIN fasta_reads on fasta_reads.sid=replnk.sid
#   WHERE map.bbb = 1 ORDER BY reftrid,origintrid};

if ( $backend eq "sqlite" ) {

    # Works close enough to the MySQL function for our purposes.
    $dbh->sqlite_create_function(
        'SUBSTRING_INDEX',
        3,
        sub {
            my ( $field, $delim, $idx ) = @_;
            my @a = split /$delim/, $field;
            return $a[$idx];
        }
    );
}

# Retrieve from the database the list of all reads mapping to a reftr called as a VNTR with at least one spanning read
my $vntr_reads_mapped_query = q{
    SELECT DISTINCT map.refid AS reftrid, SUBSTRING_INDEX(fasta_reads.head, "_", -1) AS origintrid
    FROM map INNER JOIN replnk ON replnk.rid=map.readid
        INNER JOIN fasta_reads ON fasta_reads.sid=replnk.sid
        INNER JOIN fasta_ref_reps ON map.refid = fasta_ref_reps.rid
    WHERE map.bbb = 1 AND support_vntr_span1 > 0 ORDER BY reftrid,origintrid;};
my $get_vntr_mapped_reads_sth = $dbh->prepare($vntr_reads_mapped_query);
$get_vntr_mapped_reads_sth->execute
    or die "Error executing query for all mapped reads: "
    . $get_vntr_mapped_reads_sth->errstr();
my ( $vntr, $mapped_to, %discard );
$get_vntr_mapped_reads_sth->bind_columns( \( $vntr, $mapped_to ) );

# Iterate through this list and remove all VNTRs and other mapped loci
while ( $get_vntr_mapped_reads_sth->fetch ) {
    $discard{$vntr}      = 1;
    $discard{$mapped_to} = 1;
}

for my $d ( sort keys %discard ) {
    say $d;
}

### For testing
# open my $dfile, "<", "trs.discard";

# while (my $line = <$dfile>) {
#     chomp $line;
#     $discard{$line} = 1;
# }
### End testing

my @ref_files = @run_conf{qw(REFRENCE_FILE REFERENCE_SEQ REFERENCE_INDIST)};
my $tr_count;
for my $r (@ref_files) {
    $tr_count = write_ref_file($r);
}

warn "\nDone. Final set contains $tr_count reference TRs.\n";

sub write_ref_file {
    my $file = shift;
    my @lines;
    my ( $name, $path, $fileext )
        = fileparse( $file, ".leb36", ".seq", ".indist" );

    # If not an absolute path, check in install dir
    $path = "$FindBin::Bin"
        if ( ( $path eq "./" ) && !( -e "$path/${name}${fileext}" ) );
    my $fn = "$path/${name}${fileext}";
    open my $ref_in_fh, "<", "$fn"
        or die "Error opening file $fn: $!\n";
    while ( my $line = <$ref_in_fh> ) {
        chomp $line;
        my ($trid) = split /[,\s]/, $line;

        # If this is the .seq header, we keep that.
        if ( $trid eq "Repeatid" ) {
            push @lines, $line;
            next;
        }

        push @lines, $line unless ( exists $discard{ abs($trid) } );
    }
    close $ref_in_fh;

    state $tr_count = @lines;
    state $out_suffix = ($refset_suffix) ? $refset_suffix : "refset";
    my $file_out = $tr_count . "_" . $out_suffix . $fileext;

    # Print \n after each item in @lines and at the end of the print;
    local $, = "\n";
    local $\ = "\n";
    open my $ref_out_fh, ">", $file_out;
    print $ref_out_fh @lines;
    close $ref_out_fh;

    return $tr_count;
}

# my $span1_vntrs_query = q{SELECT rid FROM fasta_ref_reps WHERE support_vntr_span1 > 0 ORDER BY rid};

# TODO Replace this in perl
# awk '{if(ARGIND==1){vntrs[$1]=1}else if((ARGIND==2)&&(-$1 in vntrs)){print $0}}' vntrs_VNTRPIPE_228486_HG38_sw100_set.span1.txt 228486_sw100_set_all_reads_mapped.txt > 228486_sw100_set_vntr_reads_mapped.txt
