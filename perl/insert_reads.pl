#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Cwd;
use DBI;
use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use lib "$RealBin/local/lib/perl5";
use Try::Tiny;

use vutil
    qw(trim get_config get_dbh set_statistics gen_exec_array_cb vs_db_insert);
use Data::Dumper;

# Arguments
my $argc = @ARGV;
die "Usage: insert_reads.pl expects 4 arguments.\n"
    unless $argc >= 4;

my $curdir          = getcwd();
my $clusterfile     = $ARGV[0];
my $indexfolder     = $ARGV[1];
my $strip454        = $ARGV[2];
my $cnf             = $ARGV[3];

my %run_conf = get_config("CONFIG", $cnf);
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

my $RECORDS_PER_INFILE_INSERT = 100000;

my %RHASH    = ();
my %HEADHASH = ();
my $count;
my $negcount = 0;
my $totalReads = 0;

my $timestart = time();

# load the RHASH now with new values added
# read clusters to see what values we will store (so we don't have to store all)
print "Reading cluster file to hash clustered ids (with added rotated repeats).\n";
open my $cluster_fh, "<", $clusterfile or die $!;
$negcount = 0;
while (<$cluster_fh>) {
    my @values = split( ',', $_ );

    foreach my $val (@values) {

        my $dir = '\'';
        if ( $val =~ /\"/ ) { $dir = '"'; }
        $val =~ s/[\'\"]//g;
        $val = trimhead($val);

        if ( $val > 0 ) {
            $RHASH{$val} = $dir;
        }
        else {
            $negcount++;
        }
    }
}
close($cluster_fh);

print keys(%RHASH) . " positive entries inserted into hash."
    . " (plus $negcount neg reference ones not in hash) ("
    . ( time() - $timestart ) . " secs).\n";

# clear  tables
print "Truncating database tables.\n";
$dbh->begin_work();
$dbh->do("DELETE FROM replnk");
$dbh->do("DELETE FROM fasta_reads");
$dbh->commit();

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");
$dbh->do("PRAGMA journal_mode = TRUNCATE");

# Insert all ref TRs from index files
my $sth = $dbh->prepare(
    qq{INSERT INTO
    replnk(rid,sid,first,last,patsize,copynum,pattern,profile,profilerc,profsize)
    VALUES (?,?,?,?,?,?,?,?,?,?)}
);

$timestart = time();
print "Reading index files and storing relevant entries in database.\n";

# KA: More file greps
opendir( DIR, $indexfolder );
my @filelist   = readdir(DIR);
my @indexfiles = sort grep( /\.index\.renumbered$/, @filelist );
my @readfiles  = sort grep( /\.reads$/, @filelist );
closedir(DIR);

my $indexcount = @indexfiles;
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
my @replnk_rows;

die "read file count doesn't equal index file count: ($indexcount vs "
    . scalar @readfiles . ")\n"
    if ( $indexcount != @readfiles );

foreach my $ifile (@indexfiles) {

    my $lfile = $ifile =~ s/index/leb36/r;

    open( $fh1, "<", "$indexfolder/$ifile" ) or die $!;
    open( $fh2, "<", "$indexfolder/$lfile" ) or die $!;

    my $line1 = read_file_line($fh1);
    my $line2 = read_file_line($fh2);

    while ( $line1 && $line2 ) {
        if ( $line1
            =~ /^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d)\t(\d+)\t([A-Z]+)$/ )
        {

            if ( exists $RHASH{"$1"} ) {
                $processed++; # this counts the number of unique ids

                $id      = $1;
                $head    = $2;
                $first   = $3;
                $last    = $4;
                $copy    = $5;
                $pat     = $6;
                $pattern = $7;

                $head = trimhead($head);

                unless ( exists $HEADHASH{"$head"} ) {
                    $HEADHASH{"$head"} = $processed;
                } # this checks for two different read ids having the same head, is that possible?

                my @values = split( ' ', $line2 );
                if ( $values[0] != $id ) {
                    die "id from index file ($id) does not match id from leb36 file ($values[0])";
                } # this makes sure the ifile and lfile are line for line parallel
                  # ... will that fail if the previous check fails?

                my $profile   = $values[5];
                my $profilerc = $values[6];
                my $proflen   = length($profile) / 2;

                push @replnk_rows,
                    [
                    $id,      $HEADHASH{"$head"}, $first,
                    $last,    $pat,               $copy,
                    $pattern, $profile,           $profilerc,
                    $proflen
                    ];

                if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                    my $cb   = gen_exec_array_cb( \@replnk_rows );
                    my $rows = vs_db_insert( $dbh, $sth, $cb,
                        "Error inserting read TRs." );
                    $inserted += $rows;
                    @replnk_rows = ();
                }
            }
        }

        $line1 = read_file_line($fh1);
        $line2 = read_file_line($fh2);

    }
    close($fh1);
    close($fh2);
}

# Remaining rows
if (@replnk_rows) {
    my $cb   = gen_exec_array_cb( \@replnk_rows );
    my $rows = vs_db_insert( $dbh, $sth, $cb, "Error inserting read TRs." );
    if ($rows) {
        $inserted += $rows;
        @replnk_rows = ();
    }
    else {
        die "\nSomething went wrong inserting, but somehow wasn't caught!\n";
    }
}

if ( $inserted == keys(%RHASH) ) {
    print $inserted . " read repeat entries inserted into database ("
        . ( time() - $timestart ) . " secs).\n";
}
else {
    die "\n\nERROR: hash contains " . keys(%RHASH)
        . " entries, while index files only have $processed matching entries and only $inserted were inserted into the database. Aborting!\n";
}

# get the read files
$totalReads = 0;
$inserted   = 0;
$processed  = 0;
$timestart  = time();
print "Reading input read files and storing relevant entries in database.\n";

$sth = $dbh->prepare(
    qq{INSERT INTO fasta_reads (sid, head, dna) VALUES(?, ?, ?)})
    or die "Couldn't prepare statement: " . $dbh->errstr;

my $headstr       = "";
my $dnastr        = "";
my $qualstr       = "";
my $HEADER_SUFFIX = "";
my @fasta_reads_rows;

my $files_processed = 0;
for my $read_file (@readfiles) {
    open my $r_fh, "<", "$indexfolder/$read_file";
    while ( my $line = <$r_fh> ) {
        chomp $line;
        ( $headstr, $dnastr ) = split "\t", $line;

        # Special last line
        # KA: This is not used anymore
        if ( $headstr eq 'totalreads' ) {
            $totalReads += $dnastr;

            # Jump out of while loop
            last;
        }
        $headstr = trimhead($headstr);
        $dnastr  = trimall($dnastr);

        my $dnabak = $dnastr;
        if ( exists $HEADHASH{$headstr} ) {

            $processed++;

            # KA: Aren't the 454 tags stripped by step 8?
            if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                warn "(insert_reads.pl) Read does not start with keyseq TCAG : "
                    . $dnabak . " (" . $headstr . ")\n";
            }

            push @fasta_reads_rows,
                [ $HEADHASH{"$headstr"}, "$headstr", "$dnastr" ];

            if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                my $cb   = gen_exec_array_cb( \@fasta_reads_rows );
                my $rows = vs_db_insert( $dbh, $sth, $cb, # this method DOES NOT RETURN if it's return value
                    "Error inserting reads.");            #   evaluates to false, and prints all the problematic lines
                $inserted += $rows;
                @fasta_reads_rows = ();
            }
        }
    }

    close $r_fh;
    $files_processed++;
}

# cleanup
if (@fasta_reads_rows) {
    my $cb   = gen_exec_array_cb( \@fasta_reads_rows );
    my $rows = vs_db_insert( $dbh, $sth, $cb,
        "Error inserting reads.");
    $inserted += $rows;
    @fasta_reads_rows = ();
}

# reenable indices
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");

# disconnect
$dbh->disconnect();

set_statistics( { NUMBER_READS => $totalReads } ) if ($totalReads);

# check ... why?

# okay, so the first two checks aren't analogous and the
# first one being true doesnt ammeliorate the second one
# being false
if ( $inserted == keys(%HEADHASH) ) {
    print $inserted . " reads inserted into database ("
        . ( time() - $timestart ) . " secs).\n";
}
elsif ( $inserted != $processed ) {
    die "\nERROR: inserted into database " . $inserted
        . " entries, while input read files have $processed matching entries. Aborting!\n";
}
else {
    die "\nERROR: hash contains " . keys(%HEADHASH)
        . " entries, while input read files only have $inserted matching entries. Aborting!\n";
}

print "Processing complete (insert_reads.pl).\n";

1;

############# subroutines ##############
sub read_file_line {
    my $fh = shift;

    if ( $fh and my $line = <$fh> ) {
        chomp $line;
        return $line;
    }
    return;
}

# Remove leading ">" if any and call trim from vutil
sub trimhead {
    my $string = shift;

    # Remove leading ">", if any
    $string =~ s/^>//;
    return trim($string);
}

# Remove all whitespace in string
sub trimall {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

# Perl flipc function to reverse direction of repeats
sub flipc {
    my $string = shift;

    ( my $ret_str = $string ) =~ tr/\'"/"\'/;
}

sub dummyquals {
    my $dna = shift;
    my @arr = split( //, $dna );
    my $len = scalar @arr;
    for ( my $i = 0; $i < $len; $i++ ) {
        $arr[$i] = 'a';
    }
    return join( "", @arr );
}
