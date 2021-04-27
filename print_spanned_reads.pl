#!/usr/bin/env perl

# prints spanned reads

sub read_file_line {
    my $fh = shift;

    if ( $fh and my $line = <$fh> ) {
        chomp $line;
        return $line;
    }
    return;
}

use strict;
use warnings;
use Cwd;
use DBI;

use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics get_trunc_query);

my $strip454 = "0";

# Perl function to remove whitespace from the start and end of the string
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Perl function to remove all whitespace
sub trimall($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

# Perl flipc function to reverse direction of repeats
sub flipc($) {
    my $string = shift;

    $string =~ s/\'/`/g;
    $string =~ s/"/\'/g;
    $string =~ s/`/"/g;

    return $string;
}

sub dummyquals($) {
    my $dna = shift;
    my @arr = split( //, $dna );
    my $len = scalar @arr;
    for ( my $i = 0; $i < $len; $i++ ) {
        $arr[$i] = 'a';
    }
    return join( "", @arr );
}

my $argc = @ARGV;

if ( $argc < 4 ) {
    die
        "Usage: print_spanned_reads.pl dbsuffix run_dir indexfolder fastafolder\n";
}

my $curdir      = getcwd;
my $DBSUFFIX    = $ARGV[0];
my $run_dir     = $ARGV[1];
my $indexfolder = $ARGV[2];
my $fastafolder = $ARGV[3];

my %run_conf = get_config( $DBSUFFIX, $run_dir );

my $totalReads = 0;

my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";
my %HEADHASH = ();

print STDERR
    "\nreading index file and storing relevant entries in database..."
    . "\n\n";
opendir( DIR, $indexfolder );
my @indexfiles = sort grep( /\.(?:index\.renumbered)$/, readdir(DIR) );
closedir(DIR);

my $indexcount = @indexfiles;
my $i          = 0;
my $inserted   = 0;
my $processed  = 0;
my $id;
my $head;
my $first;
my $last;
my $copy;
my $pat;
my $pattern;
my $fh1;
my $fh2;

foreach my $ifile (@indexfiles) {

    print STDERR "\n" . $ifile;

    open( $fh1, "<$indexfolder/$ifile" ) or die $!;
    $i = 0;

    my $line1 = read_file_line($fh1);

    while ($line1) {
        if ( $line1
            =~ /^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d)\t(\d+)\t([A-Z]+)$/ )
        {

            $id      = $1;
            $head    = $2;
            $first   = 0;
            $last    = 0;
            $copy    = 0.0;
            $pat     = 0;
            $pattern = "";

            if ( exists $HEADHASH{"$head"} ) {
            }
            else {
          #$sth0->execute("$head") or die "Cannot execute: " . $sth->errstr();
                $HEADHASH{"$head"} = $processed;
            }

            $line1 = read_file_line($fh1);
        }
    }

    close($fh1);

}

# get the fasta/fastq files
$totalReads = 0;
$inserted   = 0;
$processed  = 0;
print STDERR
    "\nreading gzipped fasta files and storing relevant entries in database..."
    . "\n\n";
opendir( DIR, $fastafolder );

# the only extensions are .tgz, .tar.gz, and .gz
#my @tarballs = grep(/\.(?:tgz|gz)$/, readdir(DIR));

my @tarballs;
@tarballs = sort grep( /^fasta.*\.(?:tgz|gz)$/, readdir(DIR) );
if ( @tarballs <= 0 ) {
    rewinddir(DIR);
    @tarballs = sort grep( /^fastq.*\.(?:tgz|gz)$/, readdir(DIR) );
}

if ( @tarballs <= 0 ) {
    die "\n\nNo fasta or fastq files found. Aborting!";
}

closedir(DIR);
my $tarball_count = @tarballs;

my $headstr       = "";
my $dnastr        = "";
my $qualstr       = "";
my $HEADER_SUFFIX = "";

FILE:
foreach (@tarballs) {
    print STDERR $_ . "";

    if ( $_ =~ /\.(?:tgz|tar\.gz)$/ ) {
        open( FASTA_IN, "tar xzfmoO '$fastafolder/$_' |" );
    }
    elsif ( $_ =~ /\.gz$/ ) {
        open( FASTA_IN, "gunzip -c '$fastafolder/$_' |" );
    }
    else {
        warn "File $_ has wrong extension. Skipping this file\n";
        next FILE;
    }

    # if this is a paired file, add .1 and .2 to all headers
    if    ( $_ =~ /^s_\d+_1/ ) { $HEADER_SUFFIX = ".1"; }
    elsif ( $_ =~ /^s_\d+_2/ ) { $HEADER_SUFFIX = ".2"; }
    else                       { $HEADER_SUFFIX = ""; }

    my $first_line = <FASTA_IN>;

    # fasta
    if ( $first_line =~ /^>(.*)(?:\n|\r)/ ) {

        #print $first_line;
        $headstr = $1;
        $dnastr  = "";
        while (<FASTA_IN>) {

            #print $_
            #end of sequence, insert and reset vars
            if ( $_ =~ /^>(.*)(?:\n|\r)/ ) {
                $headstr = trim($headstr) . $HEADER_SUFFIX;
                chomp($dnastr);
                $dnastr = trimall($dnastr);
                my $dnabak = $dnastr;
                $totalReads++;
                if ( exists $HEADHASH{"$headstr"} ) {

                    $processed++;

                    if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                        warn "Read does not start with keyseq TCAG : "
                            . $dnabak . " ("
                            . $headstr . ")\n";
                    }

                    print ">$headstr", "\n", trimall("$dnastr"), "\n";

                }

                #print "\nhead: "."`$headstr`";
                #print "\ndna: ".$dnastr;
                #last FILE;
                $headstr = $1;
                $dnastr  = "";
            }
            else {
                $dnastr .= $_;
            }
        }
        close FASTA_IN;

        # fastq
    }
    elsif ( $first_line =~ /^\@(.*)(?:\n|\r)/ ) {

        my $qmode = 0;

        #print $first_line;
        $headstr = $1;
        $dnastr  = "";
        $qualstr = "";
        while (<FASTA_IN>) {

            #print $_
            #end of sequence, insert and reset vars
            if ( $qmode == 2 && $_ =~ /^\@(.*)(?:\n|\r)/ ) {
                $qmode   = 0;
                $headstr = trim($headstr) . $HEADER_SUFFIX;
                $dnastr  = trimall($dnastr);
                my $dnabak = $dnastr;
                $qualstr = trimall($qualstr);
                $totalReads++;
                if ( exists $HEADHASH{"$headstr"} ) {

                    my $stripped = 1;

                    $processed++;

                    if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                        $stripped = 0;
                        warn "Read does not start with keyseq TCAG : "
                            . $dnabak . " ("
                            . $headstr . ")\n";
                    }

                    if ( $strip454 eq "1" && $stripped ) {

                        print ">$headstr", "\n", trimall("$dnastr"), "\n";

                    }
                    else {

                        print ">$headstr", "\n", trimall("$dnastr"), "\n";

                    }

                    $qmode = 0;
                }

                $headstr = $1;
                $dnastr  = "";
                $qualstr = "";

                # start reading quals
            }
            elsif ( $qmode == 0 && $_ =~ /^\+/ ) {
                $qmode = 1;
            }
            elsif ( $qmode != 2 ) {
                if ( $qmode == 0 ) {
                    $dnastr .= $_;
                }
                else {
                    $qualstr .= $_;
                    if ( length($qualstr) >= length($dnastr) ) { $qmode = 2; }
                }
            }
        }
        close FASTA_IN;

    }
    else {

        warn "File $_ is not a FASTA file. Skipping this file\n";
        next FILE;

    }

    # one more sequence is in the buffer, finish it off
    $headstr = trim($headstr) . $HEADER_SUFFIX;
    chomp($dnastr);
    $dnastr = trimall($dnastr);
    my $dnabak = $dnastr;
    $totalReads++;
    if ( exists $HEADHASH{"$headstr"} ) {

        my $stripped = 1;

        if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
            $stripped = 0;
            warn "Read does not start with keyseq TCAG : "
                . $dnabak . " ("
                . $headstr . ")\n";
        }

        $processed++;

        if ( $first_line =~ /^>(.*)(?:\n|\r)/ ) {
            print ">$headstr", "\n", trimall("$dnastr"), "\n";
        }
        elsif ( $first_line =~ /^\@(.*)(?:\n|\r)/ ) {
            $qualstr = trimall($qualstr);
            if ( $strip454 eq "1" && $stripped ) {
                print ">$headstr", "\n", trimall("$dnastr"), "\n";
            }
            else {
                print ">$headstr", "\n", trimall("$dnastr"), "\n";
            }
        }

    }

    print STDERR " (processed: $processed)\n";

}

# disconnect
$dbh->disconnect();

1;

