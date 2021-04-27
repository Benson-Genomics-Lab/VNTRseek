#!/usr/bin/env perl

# command line usage example:
#  ./hist_flank_filter.pl indexfolder fastafolder rotatedfolder strip_454_keytags 
#

sub read_file_line {
  my $fh = shift;

  if ($fh and my $line = <$fh>) {
    chomp $line;
    return $line;
  }
  return;
}

sub maximum ($$) { $_[$_[0] < $_[1]] }
sub minimum ($$) { $_[$_[0] > $_[1]] }

use strict;
use warnings;
use Cwd;

my $sec; 
my $min;
my $hour;
my $mday;
my $mon;
my $year;
my $wday;
my $yday;
my $isdst;

my $timestart;

my $FILTERED_BY_PATSIZE = 0;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub dummyquals($) 
{
        my $dna = shift;
	my @arr = split(//,$dna);
        my $len = scalar @arr;
        for (my $i=0; $i<$len; $i++) {
           $arr[$i] = 'a';
        }
        return join("",@arr);
}


($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf STDERR "\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

my $argc = @ARGV;


if ($argc<4) { die "Usage: hist_flank_filter.pl indexfolder fastafolder rotatedfolder strip_454_keytags\n"; }

my $curdir =  getcwd;
my $indexfolder = $ARGV[0];
my $fastafolder = $ARGV[1];
my $rotatedfolder =  $ARGV[2];
my $strip454 = $ARGV[3];

my $totalReads = 0;

my %SHASH = ();
my @HIST = ();
my %HEADHASH = ();



$timestart = time();
print STDERR "\nreading indexhist files ..."."\n\n";
opendir(DIR, $indexfolder);
my @indexfiles = grep(/\.(?:indexhist)$/, readdir(DIR));
closedir(DIR);

my $indexcount = @indexfiles;
my $i = 0;
my $j = 0;
my $inserted = 0;
my $updated = 0;
my $id; my $head; my $first; my $last; my $copy; my $pat; my $pattern;
my $fh1; my $fh2;

foreach my $ifile (@indexfiles) {

 open ($fh1, "<$indexfolder/$ifile") or die $!;
 $i=0;
 $j=0;

 print STDERR "\n".$ifile."..."; 

 my $line1 = read_file_line($fh1);

 while ($line1) {
  if ($line1 =~  /^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d)\t(\d+)\t([A-Z]+)$/ ) {

   $j++;

   { 
    $inserted++;

    $id = $1;
    $head = $2;
    $first = 0; $last = 0; $copy = 0.0; $pat = 0;    
    $pattern = "";

     {
     $first = $3;
     $last = $4;
     $copy = $5;
     $pat = $6;
     $pattern = $7;
     $i++;

     $head = trim($head);

     if (exists $HEADHASH{"$head"}) {
     } else {
      #$sth0->execute("$head") or die "Cannot execute: " . $sth->errstr();    
      $HEADHASH{"$head"} = -1;     
     }

     #print $id." ".$head." ".$first." ".$last." ".$copy." ".$pat." "." \n";

    }

   }

   #if ($i>10) { last; }
  }
  #print $_;
  $line1 = read_file_line($fh1);

 }
 close($fh1);
 print STDERR $i . " / " . $j;

}



# get the fasta/fastq files
 $totalReads=0;
 $inserted=0;
 $updated=0;
 $timestart = time();
 print STDERR "\nreading gzipped fasta files to get dna length..."."\n\n";
 opendir(DIR, $fastafolder);
 # the only extensions are .tgz, .tar.gz, and .gz
 #my @tarballs = grep(/\.(?:tgz|gz)$/, readdir(DIR));

 my @tarballs;
 @tarballs = grep(/^fasta.*\.(?:tgz|gz)$/, readdir(DIR));
 if (@tarballs <= 0) {
     rewinddir(DIR);
     @tarballs = grep(/^fastq.*\.(?:tgz|gz)$/, readdir(DIR));
 }

 if (@tarballs <= 0) {
     die "\n\nNo fasta or fastq files found. Aborting!";
 }




 closedir(DIR);
 my $tarball_count = @tarballs; 



my $headstr = "";
my $dnastr = "";
my $qualstr = "";
my $count;
my $HEADER_SUFFIX = "";

FILE:
 foreach (@tarballs) {
   print STDERR $_."";

   if ($_ =~ /\.(?:tgz|tar\.gz)$/) {
       open (FASTA_IN, "tar xzfmoO '$fastafolder/$_' |");
   } elsif ($_ =~ /\.gz$/) {
       open (FASTA_IN, "gunzip -c '$fastafolder/$_' |");
   } else {
       warn "File $_ has wrong extension. Skipping this file\n";
       next FILE;
   }

   # if this is a paired file, add .1 and .2 to all headers
   if ($_ =~ /^s_\d+_1/) { $HEADER_SUFFIX = ".1"; }
   elsif ($_ =~ /^s_\d+_2/) { $HEADER_SUFFIX = ".2"; }
   else { $HEADER_SUFFIX = ""; }

    my $first_line = <FASTA_IN>;

    # fasta
    if ($first_line =~ /^>(.*)(?:\n|\r)/) {

                                #print $first_line;
				$headstr = $1;
				$dnastr="";
                                while (<FASTA_IN>) { 
				 #print $_ 
				 #end of sequence, insert and reset vars
				 if ($_ =~ /^>(.*)(?:\n|\r)/) {
				  $headstr = trim($headstr).$HEADER_SUFFIX;
                                  chomp($dnastr);
				  $dnastr = trim($dnastr);
				  my $dnabak = $dnastr;
				  $totalReads++;
				  if (exists $HEADHASH{"$headstr"}) {


                                      if ($strip454 eq "1" && $dnastr !~ s/^TCAG//i) {
                                                        warn "Read does not start with keyseq TCAG : ".$dnabak." (".$headstr.")\n";
				      }

			              $inserted++;
				      $HEADHASH{"$headstr"} = length("$dnastr");

		  		      if (0 == $HEADHASH{"$headstr"}) { warn "\n length of '$headstr' is ". $HEADHASH{"$headstr"} . "! (warning 1)\n"; }


				  }
				  #print "\nhead: "."`$headstr`";
				  #print "\ndna: ".$dnastr;
				  #last FILE;
				  $headstr = $1;
				  $dnastr="";
				 } else {
  				  $dnastr.=$_;
				 }
				}
                                close FASTA_IN;

       # fastq
       } elsif ($first_line =~ /^\@(.*)(?:\n|\r)/) {

				my $qmode=0;

                                #print $first_line;
				$headstr = $1;
				$dnastr="";
				$qualstr="";
                                while (<FASTA_IN>) { 
				 #print $_ 
				 #end of sequence, insert and reset vars
				 if ($qmode==2 && $_ =~ /^\@(.*)(?:\n|\r)/) {
				  $qmode=0;
				  $headstr = trim($headstr).$HEADER_SUFFIX;
				  $dnastr = trim($dnastr);
				  my $dnabak = $dnastr;
				  $qualstr = trim($qualstr);
				  $totalReads++;
				  if (exists $HEADHASH{"$headstr"}) {

				      my $stripped = 1;

                                      if ($strip454 eq "1" && $dnastr !~ s/^TCAG//i) {
							$stripped = 0;
                                                        warn "Read does not start with keyseq TCAG : ".$dnabak." (".$headstr.")\n";
				      } 
 			              $inserted++;

				      $HEADHASH{"$headstr"} = length("$dnastr");

		 		      if (0 == $HEADHASH{"$headstr"}) { warn "\n length of '$headstr' is ". $HEADHASH{"$headstr"} . "! (warning 2)\n"; }


				      $qmode=0;
				  } 

				  $headstr = $1;
				  $dnastr="";
				  $qualstr="";

				 # start reading quals
				 } elsif ($qmode==0 && $_ =~ /^\+/) {
				    $qmode=1;
				 } elsif ($qmode != 2) {
   				     if ($qmode == 0)  {
  				       $dnastr.=$_;
				     } else {
   				       $qualstr.=$_;
				       if (length($qualstr) >= length($dnastr)) { $qmode = 2; }
				     }
				 }
				}
                                close FASTA_IN;


       } else {
           
		warn "File $_ is not a FASTA file. Skipping this file\n";
                next FILE;

	}


       # one more sequence is in the buffer, finish it off
       $headstr = trim($headstr).$HEADER_SUFFIX;
       chomp($dnastr);
       $dnastr = trim($dnastr);
       my $dnabak = $dnastr;
       $totalReads++;
       if (exists $HEADHASH{"$headstr"}) {

 	        my $stripped = 1;

                if ($strip454 eq "1" && $dnastr !~ s/^TCAG//i) {
					$stripped = 0;
                                        warn "Read does not start with keyseq TCAG : ".$dnabak." (".$headstr.")\n";
                }

                $inserted++;


                if ($first_line =~ /^>(.*)(?:\n|\r)/) {
 		        $HEADHASH{"$headstr"} = length("$dnastr");
                } elsif ($first_line =~ /^\@(.*)(?:\n|\r)/) {
	 	        $HEADHASH{"$headstr"} = length("$dnastr");
                }

		if (0 == $HEADHASH{"$headstr"}) { warn "\n length of '$headstr' is ". $HEADHASH{"$headstr"} . "! (warning 3)\n"; }


       }



 print STDERR " (reads with TRs: $inserted)\n";

 }


print STDERR "\ncounting filtered and producing histogram..."."\n\n";
foreach my $ifile (@indexfiles) {

 open ($fh1, "<$indexfolder/$ifile") or die $!;
 $i=0;
 $j=0;

 print STDERR "\n".$ifile."...";

 my $line1 = read_file_line($fh1);

 while ($line1) {
  if ($line1 =~  /^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d)\t(\d+)\t([A-Z]+)$/ ) {

   $j++;

   {
    $inserted++;

    $id = $1;
    $head = $2;
    $first = 0; $last = 0; $copy = 0.0; $pat = 0;
    $pattern = "";

     {
     $first = $3;
     $last = $4;
     $copy = $5;
     $pat = $6;
     $pattern = $7;
     $i++;

     $head = trim($head);

       if (exists $HEADHASH{"$head"} && -1 != $HEADHASH{"$head"}) {
      	 my $mflank = minimum( ($first-1), ($HEADHASH{"$head"} - $last));

         if ($pat>=7 && $mflank<=20) {
	   $HIST[$mflank]++;
         }

	 if ($pat<7) {
	   $FILTERED_BY_PATSIZE++;
	 }
       }
     
     #}

     #print $id." ".$head." ".$first." ".$last." ".$copy." ".$pat." "." \n";

    }

   }

   #if ($i>10) { last; }
  }
  #print $_;
  $line1 = read_file_line($fh1);

 }
 close($fh1);
 print STDERR $i . " / " . $j;

}


print STDERR "\n\nFiltered by patsize: $FILTERED_BY_PATSIZE\n\n";
print STDERR "\n\nFiltered by flanks:\n\n";

for ($i=0; $i<=20; $i++) { 
 my $sum=0;
 for ($j=0; $j<=$i; $j++) { 
    $sum += $HIST[$j];
  }
  print STDERR $i . ". " . $sum . "\n";
}

print STDERR "\n\nProcessing complete (hist_flank_filter.pl).\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf STDERR "\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;




1;




