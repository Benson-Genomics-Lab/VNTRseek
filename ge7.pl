#!/usr/bin/env perl

use strict;
use warnings;

my $TRCOUNTge7 = 0;
my $READCOUNTge7 = 0;
my $READCOUNT = 0;

my %HASH = ();
my %HASHALL = ();

while (<>) {
 my @l = split('\t',$_);

 if ($l[-2]>=7) { 
   $TRCOUNTge7++;
   my $header = $l[1];
   if (!$HASH{$header}) {
     #print $header;
     $HASH{$header} = 1;
     $READCOUNTge7++;
   }
 }


 {
   my $header = $l[1];
   if (!$HASHALL{$header}) {
     #print $header;
     $HASHALL{$header} = 1;
     $READCOUNT++;
   }
 }


}

print "$TRCOUNTge7 $READCOUNTge7 $READCOUNT";
