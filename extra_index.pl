#!/usr/bin/env perl

# prints sequences for BBB for pcr_dup step

use strict;
use warnings;
use Cwd;
use POSIX qw(strftime);
use DBI;
use File::Basename;

use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics get_trunc_query);

my $FASTA = 1;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $argc = @ARGV;

if ( $argc < 3 ) { die "Usage: extra_index.pl expects 3 arguments!\n"; }

my $curdir = getcwd;

my $folder   = $ARGV[0];
my $DBSUFFIX = $ARGV[1];
my $run_dir  = $ARGV[2];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $DBSUFFIX, $run_dir );
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

# create folder
my $exstring = "rm -rf $folder";
system($exstring);
mkdir($folder);

my $sth;
my $query;
my $result;
my $num;
my $i;

my $sql_clause = q{
  FROM map INNER JOIN
    rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN
    rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid INNER JOIN
    replnk ON replnk.rid=map.readid INNER JOIN
    fasta_reads on fasta_reads.sid=replnk.sid
  ORDER BY map.refid,map.readid};
($num) = $dbh->selectrow_arrayref( q{SELECT COUNT(*) } . $sql_clause )
    or die "Couldn't execute statement: " . $dbh->errstr;
$sth
    = $dbh->prepare(
    q{SELECT map.refid, map.readid, replnk.sid, replnk.first, replnk.last, replnk.copynum, replnk.patsize, replnk.pattern,fasta_reads.dna}
        . $sql_clause )
    or die "Couldn't prepare statement: " . $dbh->errstr;

$sth->execute() or die "Couldn't execute: " . $sth->errstr;

print "\n best best best records: $num\n";

my $oldref  = -1;
my $oldread = -1;
$i = 0;
my $nrefs = 0;

my $fh;
while ( my @data = $sth->fetchrow_array() ) {
    if ( $data[0] != $oldref ) {
        if ( $i != 0 ) { close($fh); }
        $nrefs++;
        print "\n$nrefs";
        open $fh, ">$folder/$data[0].seq" or die $!;
    }
    print $fh "$data[1] $data[2] $data[3] $data[4]";
    printf $fh " %.2lf", $data[5];
    $data[8] =~ s/\s+//g;
    print $fh " $data[6] $data[7] $data[8]\n";
    $oldref  = $data[0];
    $oldread = $data[1];
    $i++;
}

close($fh) if ($fh);

$sth->finish;

print "\n\nProcessing complete (extra_index.pl), $nrefs files created.\n";

warn strftime( "\n\nend: %F %T\n\n", localtime );
#
$dbh->disconnect();

1;

