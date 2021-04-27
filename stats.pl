#!/usr/bin/env perl

use warnings;

use DBI;
use List::Util qw[min max];

use FileHandle;
use Getopt::Std;
use File::Copy;

#use 5.012;

my $DBNAME = "V";
my $LOGIN = "";
my $PASS = "";
my $HOST = "localhost";
my $sth;
my $sth2;

my $dbh = DBI->connect("DBI:mysql:$DBNAME;host=$HOST", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr";


my $NLINES = 0;
my %HEADHASH = ();
my %CLUSTEREDHASH = ();

my $RPINTERSECT = 0;
my $NEITHERONLY = 0;
my $NEITHERONLY_TRS2ORMORE = 0;
my $PROFONLY = 0;
my $PROFONLY_TRS2ORMORE = 0;
my $FLANKONLY = 0;
my $FLANKONLY_TRS2ORMORE = 0;
my $PTIE = 0;
my $FTIE = 0;
 
 open (FILE,"spanning_reads_not_mapped_trs.index") or die "Coud not open input file!\n";
 while (<FILE>) {
 chomp;

  if (m/^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d+)/) {

    $NLINES++;
   
    # unique reads
     if (exists $HEADHASH{"$2"}) {
     } else {
      $HEADHASH{"$2"}= 1;
     }

#   print "result: $1|$2";

#   $sth = $dbh->prepare("select head from fasta_reads inner join replnk ON replnk.sid=fasta_reads.sid where rid=$1") or die "Couldn't prepare statement: " . $dbh->errstr;
#    $sth->execute();
#    if ($sth->rows > 0) {
#    my @data = $sth->fetchrow_array();
#    print " \"".$data[0]."\"";
#    }
#    $sth->finish();
   
#   print "\n";

#   exit(0);
  }

 }
 close (FILE);

 # clustered?
 $sth = $dbh->prepare("select rid from fasta_reads inner join replnk ON replnk.sid=fasta_reads.sid where head=?") or die "Couldn't prepare statement: " . $dbh->errstr;
 foreach my $key(keys(%HEADHASH)) {
    $sth->execute($key);
    if ($sth->rows > 0) {
    my @data = $sth->fetchrow_array();
      $CLUSTEREDHASH{"$key"}=1;
      # print $key." ".$data[0]."\n";  
    }
 }
 $sth->finish();

 # intestect and ties (could produce multipele records but the only affected thing would be random assignment to profonly or flankonly)
 $sth = $dbh->prepare("select replnk.sid,rank.ties,rankflank.ties,rank.readid,rankflank.readid,rank.refid,rankflank.refid,map.readid,map.refid from map LEFT JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid LEFT JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid INNER JOIN replnk ON replnk.rid=map.readid INNER JOIN fasta_reads on fasta_reads.sid=replnk.sid WHERE head=?") or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth2 = $dbh->prepare("select count(*) from replnk where sid=?");
 my $count=0;
 foreach my $key(keys(%CLUSTEREDHASH)) {
    $sth->execute($key);
    my $sid = 0;
    my $counttrs = 0;
    my $isptie = 0;
    my $isftie = 0;
    my $isinrank = 0; my $rankid = 0; my $frankid = 0; my $rankidR = 0; my $frankidR = 0; my $mapTR = 0; $mapTRref = 0;
    my $isinrankflank = 0;
    my $num = 0;
    my $rows = $sth->rows;
    while ($num < $rows) {
        
      $mapTRref = 0; $mapTR = 0; $sid = 0; $isptie = 0; $isftie = 0; $isinrank = 0; $rankid = 0; $rankidR = 0; $isinrankflank = 0; $frankid = 0; $frankidR = 0;
      my @data = $sth->fetchrow_array();
      if (defined $data[0] && $data[0]>0) { $sid=$data[0]; }
      if (defined $data[1] && $data[1]>0) { $isptie=1; }
      if (defined $data[2] && $data[2]>0) { $isftie=1; }
      if (defined $data[3] && $data[3]>0) { $isinrank=1; $rankid=$data[3]; $rankidR=$data[5]; }
      if (defined $data[4] && $data[4]>0) { $isinrankflank=1; $frankid=$data[4]; $frankidR=$data[6]; }
      if (defined $data[7] && $data[7]>0) { $mapTR=$data[7]; }
      if (defined $data[8] && $data[8]>0) { $mapTRref=$data[8]; }

      # base stats on first left left join entry for now
      if ($num == 0) {

        $count++;

        $sth2->execute($sid);
        if ($sth2->rows > 0) {
          my @data = $sth2->fetchrow_array();
          $counttrs = $data[0];
        }

        if ($isinrank && $isinrankflank) { $RPINTERSECT++; } 
          else 
        { 
           if ($isinrank || $isinrankflank) { 
                if ($isinrank) { $PROFONLY++;  }
                if ($isinrankflank) { $FLANKONLY++;  }
                if ($counttrs>=2 && $isinrank) { $PROFONLY_TRS2ORMORE++;} 
                if ($counttrs>=2 && $isinrankflank) { $FLANKONLY_TRS2ORMORE++;} 
          } else { 
               $NEITHERONLY++; if ($counttrs>=2) { $NEITHERONLY_TRS2ORMORE++;} 
          } 
       }
       if ($isinrank && $isinrankflank && $isptie) { $PTIE++; }
       if ($isinrank && $isinrankflank && $isftie) { $FTIE++; }


      } # end of checking for first

     { print "$count.$num ($key) (isinrank-isinrankflank-ptie-ftie) = $isinrank-$isinrankflank-$isptie-$isftie trs/read=$counttrs\n\tMAP=$mapTR->$mapTRref P=($rankid -> $rankidR) F=($frankid -> $frankidR) \n";}

      $num++;

    } # end of check sid key

    print "\n";

#    if ($count>=200) { last; }
 }
 $sth->finish();
 $sth2->finish();

 # results
 print "\n\n";
 print "INPUT LINES: $NLINES\n";
 print "UNIQUE READS: " .  keys(%HEADHASH) . " \n";
 print "READS CLUSTERED: ". keys(%CLUSTEREDHASH) .  " \n";
 print "READS PROF/FLANK INTERSECT (ties?): ". $RPINTERSECT .  " \n";
 print "\tprofile tie: ". $PTIE .  " \n";
 print "\tflank tie: ". $FTIE .  " \n";
 print "READS PROF ONLY : ". $PROFONLY . ", 2 or more trs: ". $PROFONLY_TRS2ORMORE. " \n";
 print "READS FLANK ONLY : ". $FLANKONLY . ", 2 or more trs: ". $FLANKONLY_TRS2ORMORE. " \n";
 print "READS NEITHER PROF OR FLANK: ". $NEITHERONLY . ", 2 or more trs: ". $NEITHERONLY_TRS2ORMORE.  " \n";


 exit;
