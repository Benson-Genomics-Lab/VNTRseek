#!/usr/bin/env perl

# command line usage example:
#  ./checkleb36.pl inputfolder

use strict;
use warnings;
use Cwd;

die "Useage: checkl3b36.pl expects 1 argument.\n"
    unless @ARGV;

my $curdir =  getcwd();
my $tgz_dir = $ARGV[0];

# KA: Bad file grep for list.
# get a list of input files
opendir(DIR, $tgz_dir);
# the only extensions are .leb36
my @tarballs = grep(/\.(?:leb36)$/, readdir(DIR));
closedir(DIR);

my $tarball_count = @tarballs;
print "$tarball_count supported files found in $tgz_dir\n";
die "Exiting\n" if $tarball_count == 0;

# enter dir
chdir($tgz_dir);

print "Checking read leb36 files.\n";
my %uhash = ();
my $fh;
foreach my $ifile (@tarballs) {
    open($fh, "<$ifile") or die $!;
    while (<$fh>) {
        if (/^(\d+)/) {
            if (exists $uhash{$1}) {
                die "Non-unique id ($1) detected in reads."
                    . " Were steps 2 and 3 executed?\n";
            }
            $uhash{$1} = 1;
        }
    }
    close $fh;
}
