#!/usr/bin/env perl

use strict;
use warnings;

use DBI;
use List::Util qw[min max];

use FileHandle;
use Getopt::Std;
use File::Copy;


# Prints references that are not clustered


if (@ARGV<3) {
 print STDERR "\n\nprint_unclustered.pl: Enter database name, login and pass!\n";
 exit 1;
}

my $DBNAME = $ARGV[0];
my $LOGIN = $ARGV[1];
my $PASS = $ARGV[2];
my $i;
my $num;
my $sth;

my $dbh = DBI->connect("DBI:mysql:$DBNAME", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr";




# populate a hash of clustered
my %MYHASH = ();
$sth = $dbh->prepare('SELECT repeatid from clusterlnk WHERE repeatid<0;;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$num = $sth->rows;
for ($i=0; $i<$num; $i++) {
   my @data = $sth->fetchrow_array();
   my $val = -1 * $data[0];
   $MYHASH{$val}=1;
   #print "$val\n";
}
$sth->finish();



# get all not in hash
$sth = $dbh->prepare('SELECT rid,sequence,flankleft,flankright from fasta_ref_reps;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$num = $sth->rows;
for ($i=0; $i<$num; $i++) {
   my @data = $sth->fetchrow_array();
   my $val = int($data[0]);
   if (exists $MYHASH{$val}) {
     #print "removing $val!\n";
   } else {
     print ">".$data[0]."\n". lc(substr($data[2],-20)). uc($data[1]). lc(substr($data[3],0,20)). "\n";
   }
}
$sth->finish();




