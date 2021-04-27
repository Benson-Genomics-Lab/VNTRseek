#!/usr/bin/env perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use List::Util qw[min max];

# this is similar to run_variablility.pl, except it's more liberal at picking targets and updates
# different fields
# command line usage example:
#  ./run_assemblyreq.pl inputfile mapdir dbname dblogin dbpass
# where inputfile is the main cluster file
#

use strict;
use warnings;
use Cwd;
use DBI;

use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib"; 
require "vutil.pm";

use vutil ('get_credentials');
use vutil ('write_mysql');
use vutil ('stats_set');

my $updatedClustersCount = 0;


my $sec; 
my $min;
my $hour;
my $mday;
my $mon;
my $year;
my $wday;
my $yday;
my $isdst;

sub nowhitespace($)
{
	my $string = shift;
	$string =~ s/\s+//g;
	return $string;
}


($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);

my $curdir =  getcwd;

my $argc = @ARGV;
if ($argc<6) { die "Usage: run_assemblyreq.pl inputfile mapdir dbname msdir minflank tempdir\n"; }


my $inputfile = $ARGV[0];
my $mapdir = $ARGV[1];
my $DBNAME = $ARGV[2];
my $MSDIR = $ARGV[3];
my $MIN_FLANK_REQUIRED = $ARGV[4];
my $TEMPDIR = $ARGV[5];

# set these mysql credentials in vs.cnf (in installation directory)
my ($LOGIN,$PASS,$HOST) = get_credentials($MSDIR);


####################################
sub SetStatistics {

  my $argc = @_;
  if ($argc <2) { die "stats_set: expects 2 parameters, passed $argc !\n"; }

  my $NAME = $_[0];
  my $VALUE = $_[1];

  #print "$DBNAME,$LOGIN,$PASS,$NAME,$VALUE\n";
  return stats_set($DBNAME,$LOGIN,$PASS,$HOST,$NAME,$VALUE);
}
####################################



my $clusters_processed = 0;


my $dbh = DBI->connect("DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST", "$LOGIN", "$PASS"
	           ) || die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth1;
my $sth2;
my $sth3;
my $sth5;
my $sth7;
my $sth8;

# change settings to speedup updates and inserts
$sth = $dbh->prepare('SET AUTOCOMMIT = 0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET FOREIGN_KEY_CHECKS = 0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET UNIQUE_CHECKS = 0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;


# update reserved field on entire table
$sth = $dbh->prepare('UPDATE clusterlnk SET reserved2=0')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;
$sth = $dbh->prepare('UPDATE map SET reserved2=0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;


 # create temp table for updates
 $sth = $dbh->prepare('CREATE TEMPORARY TABLE mapr (`refid` INT(11) NOT NULL, `readid` INT(11) NOT NULL, PRIMARY KEY (refid,readid)) ENGINE=INNODB;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE mapr DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;


# prepare statments
$sth = $dbh->prepare('SELECT rid,flankleft,sequence,flankright,pattern,copynum,(lastindex-firstindex+1) as arlen FROM fasta_ref_reps WHERE rid = ?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth1 = $dbh->prepare('SELECT rid, dna, first, last, pattern, copynum FROM fasta_reads INNER JOIN replnk ON fasta_reads.sid=replnk.sid WHERE rid = ?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth2 = $dbh->prepare('UPDATE clusters SET assemblyreq=? WHERE cid=?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth3 = $dbh->prepare('UPDATE clusterlnk SET reserved2=1 WHERE clusterid=? AND repeatid=?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth5 = $dbh->prepare('SELECT SUM(reserved2) FROM clusterlnk WHERE clusterid=?;')
                or die "Couldn't prepare statement: " . $dbh->errstr;

#$sth7 = $dbh->prepare('UPDATE map SET reserved2=1 WHERE refid=? AND readid=?')
$sth7 = $dbh->prepare("LOAD DATA LOCAL INFILE '$TEMPDIR/mapr_$DBNAME.txt' INTO TABLE mapr FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;

$sth8 = $dbh->prepare('SELECT map.refid, map.readid from clusterlnk INNER JOIN map ON map.refid=-clusterlnk.repeatid WHERE map.bbb=1 AND clusterlnk.clusterid=? ORDER BY map.readid,map.refid;')
                or die "Couldn't prepare statement: " . $dbh->errstr;

my $processed=0;

my $TEMPFILE;
open ($TEMPFILE, ">$TEMPDIR/mapr_$DBNAME.txt") or die $!;

# this is a function that detects if reads have (somewhat) different number of copy numbers then reference(s)
sub VNTR_YES_NO {

 my(%refs) = %{$_[0]};
 my(%reads) = %{$_[1]};
 my(%readhash) = %{$_[2]};
 my $clusterid = $_[3];

 my $varyes = 0;

 my(%REF_UPDATED) = ();

 #print "\nVNTR_YES_NO:";

 # for each read in map file
 while ( my ($key, @temp) = each(%readhash) ) {

	my $valread;
	my $valref;

   if (exists $reads{$key}) {
 
 	$valread = $reads{$key};

        #print "\n$key($valread) =>";

	# for each references listed in map file as being associated with the read
        foreach my $val (@{ $readhash{$key}  }) {

	  if (exists $refs{$val}) {
		$valref = $refs{$val};
	  } else {
		print "Undefined reference `$val`. Aborting!\n";
            	exit(1);
  	  }

  	  #print " ".$val."($valref)";

	  # print "$valRead - $valRef\n";
    	  #if (($valread >= ($valref-.5)) && ($valread <= ($valref+.5))) {
	  {
		$varyes = 1;

		#print "\t$clusterid: $refid\n";		

		if (!exists $REF_UPDATED{$val}) {

  		  $sth3->execute($clusterid,$val)             # set the reserved field on reference
	             or die "Couldn't execute statement: " . $sth3->errstr;

	 	  $REF_UPDATED{$val}=1;
		}

		#$sth7->execute(-$val,$key)             # set the reserved field on map read, (added 11/19/2010)
	        #    or die "Couldn't execute statement: " . $sth7->errstr;

		$processed++;

		print $TEMPFILE -$val,",",$key,"\n";

	        if ($processed % $RECORDS_PER_INFILE_INSERT == 0) {
		       close($TEMPFILE);
		       $sth7->execute();
		       open ($TEMPFILE, ">$TEMPDIR/mapr_$DBNAME.txt") or die $!;
     		}
		

	  }
	
	}
	#print "\n";


   } # end of exists

 } # end of while loop

 # otherwise mark as not variable
 # print "\nNOT VARIABLE\n\n";
 return $varyes;

} # end of func





open FILE, "<$inputfile" or die $!;



while (<FILE>) { 
 $clusters_processed++;

 my @values = split(',', $_);

 my %REFCOPIES = ();
 my %ASKLENGTH = ();
 my %READCOPIES = ();
 my $repeatcount = 0;
 my $refcount = 0;
 my $readcount = 0;

 my $readlen;
 my $first;
 my $last;
 my $vntrres;

 foreach my $val (@values) {

   my $dir = '\'';
   if ($val =~ m/([\'\"])/) { $dir = $1; } 

   chomp($val);
 
   $val =~ s/[\'\"]//g;

   $repeatcount++;

   # go though all refs for each read
   if ($val<=0) {

     $refcount++;	

     #$val = ($val<0)? -$val : $val; # ids are positive in database

     $sth->execute(-$val)             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

     my @data = $sth->fetchrow_array();

     if ($sth->rows == 0) {
            print STDERR "No record in database for entry `$val'. Aborting!\n\n";
 	    exit(1);
     }
   

     $ASKLENGTH{$val} = int($data[6]); # remember length of the ref array sequence
     $REFCOPIES{$val} = $data[5];
  
   } else {

     $readcount++;

     $sth1->execute($val)             # Execute the query
            or die "Couldn't execute statement: " . $sth1->errstr;

     my @data = $sth1->fetchrow_array();

     if ($sth1->rows == 0) {
            print STDERR "No record in database for entry `$val'. Aborting!\n";
 	    exit(1);
     }

     $readlen = length(nowhitespace($data[1]));

     $ASKLENGTH{$val} = $readlen; # remember length of the read

     $first = $data[2];
     $last = $data[3];
     if (($first-1) < $MIN_FLANK_REQUIRED || ($readlen - $last) < $MIN_FLANK_REQUIRED) {
        $READCOPIES{$val} = $data[5];
     }
   }

 }

 # get map data
 my %READVECTOR = ();
 my $readidold=-1;
 $sth8->execute($clusters_processed)             # Execute the query
            or die "Couldn't execute statement: " . $sth8->errstr;
 my $nmapped = $sth8->rows;
 for (my $i=0; $i<$nmapped; $i++) {
     my @data = $sth8->fetchrow_array();

     my $readid;
     my $refid;
     my $readlen;

     $refid = -1 * int($data[0]);
     $readid = int($data[1]);
    
print STDERR "$refid | $readid \n";

     if ($readid != $readidold) {
       $READVECTOR{$readid} = ();
     }     

     $readlen = $ASKLENGTH{$readid};
     #if ($ASKLENGTH{$refid} <= $readlen) {
       push( @{ $READVECTOR{$readid}  }, $refid  );
     #}

     $readidold=$readid;
 }
$sth8->finish;
 print STDERR "\nprocessed: $clusters_processed\n\n";
 

 # do for 1st 10 clusters for now
 # if ($clusters_processed >= 20) { last; }

 # insert database records (cluster table)
 #print "VYN:  ";
 #print  VNTR_YES_NO(\%REFCOPIES,\%READCOPIES,\%READVECTOR);

 my $vYes=VNTR_YES_NO(\%REFCOPIES,\%READCOPIES,\%READVECTOR,$clusters_processed);


 # get number of updated refs
 $sth5->execute($clusters_processed)             # Execute the query
            or die "Couldn't execute statement: " . $sth5->errstr;
 my @data = $sth5->fetchrow_array();
 $vntrres=0;
 if ($sth5->rows != 0) { $vntrres=$data[0]; }

 $updatedClustersCount += $vYes;

 $sth2->execute($vntrres,$clusters_processed)  # Execute the query
         or die "Couldn't execute statement: " . $sth2->errstr;



}
close(FILE);


 # finish loading the file into tempfile and switch indexes back on
 close($TEMPFILE);
 $sth7->execute();
 $sth7->finish;

 $sth7 = $dbh->prepare('ALTER TABLE mapr ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth7->execute();
 $sth7->finish;


 # update based on temp table
 my $updfromtable = 0;
 my $query = "UPDATE mapr p, map pp SET reserved2=1 WHERE pp.refid = p.refid AND pp.readid = p.readid;";
 $sth7 = $dbh->prepare($query);
 $updfromtable=$sth7->execute();
 $sth7->finish;


 # cleanup temp file
 unlink("$TEMPDIR/mapr_$DBNAME.txt");



$sth->finish;
$sth1->finish;
$sth2->finish;
$sth3->finish;
$sth5->finish;
#$sth7->finish;


 # set old settings
$sth = $dbh->prepare('SET AUTOCOMMIT = 1;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET FOREIGN_KEY_CHECKS = 1;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET UNIQUE_CHECKS = 1;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;






$dbh->disconnect();


if ($updfromtable != $processed) { die "Updated number of entries($updfromtable) not equal to the number of updated counter ($processed), aborting!"; }

print STDERR "\n\n";

print STDERR "Processing complete -- processed $clusters_processed cluster(s).\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);



1;




