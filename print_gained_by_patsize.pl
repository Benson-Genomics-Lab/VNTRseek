#!/usr/bin/env perl

use strict;
use warnings;

use DBI;
use List::Util qw[min max];

use FileHandle;
use Getopt::Std;
use File::Copy;


# Prints references that are not clustered

my $MIN_SUPPORT_REQUIRED = 2;

if (@ARGV<5) {
 print STDERR "\n\nprint_unclustered.pl: Enter database name, login and pass and MINPAT and MAXPAT!\n";
 exit 1;
}

my $DBNAME = $ARGV[0];
my $LOGIN = $ARGV[1];
my $PASS = $ARGV[2];
my $MINPAT = $ARGV[3];
my $MAXPAT = $ARGV[4];
my $i;
my $num;
my $sth;

my $dbh = DBI->connect("DBI:mysql:$DBNAME", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr";


# copies gained/lost
my $total = 0;
print "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) PatternSize, Copies Gained, Frequency, ArraySize\n";
$sth = $dbh->prepare("SELECT length(pattern) AS PSIZE, copies, count(*),(lastindex-firstindex+1) AS arraysize FROM vntr_support INNER JOIN fasta_ref_reps ON vntr_support.refid=-fasta_ref_reps.rid AND (lastindex-firstindex+1)>=$MINPAT AND (lastindex-firstindex+1)<=$MAXPAT WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED GROUP BY PSIZE ASC, copies ASC;")
                or die "Couldn't prepare statement: " . $dbh->errstr;
$num = 0;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
$num = $sth->rows;
$i=0;
while ($i<$num) {
  my @data = $sth->fetchrow_array();
  print $data[0] . "," .  $data[1] . "," .  $data[2] . "," .  $data[3] . "\n";
  $i++;
  $total += $data[2];
}
$sth->finish;
print "TOTAL: $total\n";

