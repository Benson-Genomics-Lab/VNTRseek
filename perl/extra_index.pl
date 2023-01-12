#!/usr/bin/env perl

# prints sequences for BBB for pcr_dup step

use strict;
use warnings;
use Cwd;
use DBI;
use File::Basename;

use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics);

# Arguments
my $argc = @ARGV;
die "Usage: extra_index.pl expects 2 arguments.\n"
    unless $argc >= 2;

my $curdir = getcwd();
my $folder = $ARGV[0];
my $cnf    = $ARGV[1];

# database
my %run_conf = get_config("CONFIG", $cnf);
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

my $sth;

my $sql_clause = q{
  FROM map JOIN rank using(refid, readid)
      JOIN rankflank using(refid, readid)
      JOIN replnk ON replnk.rid=map.readid
      JOIN fasta_reads on fasta_reads.sid=replnk.sid
  ORDER BY map.refid, map.readid};
my ($num) = $dbh->selectrow_array( q{SELECT COUNT(*) } . $sql_clause )
    or die "Couldn't execute statement: " . $dbh->errstr;
$sth = $dbh->prepare(
    q{SELECT map.refid, map.readid, replnk.sid, replnk.first, replnk.last,
      replnk.copynum, replnk.patsize, replnk.pattern, fasta_reads.dna}
        . $sql_clause )
    or die "Couldn't prepare statement: " . $dbh->errstr;

$sth->execute() or die "Couldn't execute: " . $sth->errstr;

print "Best best best records: $num\n";

my $oldref  = -1;
my $i = 0;
my $nrefs = 0;

my $fh;
while ( my @data = $sth->fetchrow_array() ) {
    if ( $data[0] != $oldref ) {
        if ( $i != 0 ) { close($fh); }
        $nrefs++;
        open $fh, ">$folder/$data[0].seq" or die $!;
    }
    print $fh "$data[1] $data[2] $data[3] $data[4]";
    printf $fh " %.2lf", $data[5];
    $data[8] =~ s/\s+//g;
    print $fh " $data[6] $data[7] $data[8]\n";
    $oldref  = $data[0];
    $i++;
}
close($fh) if ($fh);
$sth->finish();
$dbh->disconnect();

print "Processing complete (extra_index.pl), $nrefs files created.\n";
