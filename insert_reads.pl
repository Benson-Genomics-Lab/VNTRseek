#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Cwd;
use DBI;
use File::Basename;
use POSIX qw(strftime);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use lib "$RealBin/local/lib/perl5";
use Try::Tiny;

# use ProcInputReads
# qw(get_reader init_bam formats_regexs compressed_formats_regexs set_install_dir);
# set_install_dir("$FindBin::RealBin");

use vutil
    qw(get_config get_dbh set_statistics get_trunc_query gen_exec_array_cb vs_db_insert);
use Data::Dumper;

my $RECORDS_PER_INFILE_INSERT = 100000;

my $timestart;

my $count;

################ main ####################
warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $argc = @ARGV;

if ( $argc < 10 ) {
    die
        "Usage: insert_reads.pl  clusterfile indexfolder fastafolder rotatedfolder rotatedreffile strip_454_keytags dbname msdir tempdir ispairedreads\n";
}

my $curdir          = getcwd;
my $clusterfile     = $ARGV[0];
my $indexfolder     = $ARGV[1];
my $fastafolder     = $ARGV[2];
my $rotatedfolder   = $ARGV[3];
my $rotatedreffile  = $ARGV[4];
my $strip454        = $ARGV[5];
my $DBSUFFIX        = $ARGV[6];
my $run_dir         = $ARGV[7];
my $TEMPDIR         = $ARGV[8];
my $IS_PAIRED_READS = $ARGV[9];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $DBSUFFIX, $run_dir );

my $totalReads = 0;

# TODO for all files needing this function, maybe run get_config first
# to eliminate need for second arg
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";
my %RHASH    = ();
my %SHASH    = ();
my %HEADHASH = ();
my $negcount = 0;

$timestart = time();

system("cp $clusterfile $rotatedfolder/allwithdups.clusters");

# load the RHASH now with new values added
#
# read clusters to see what values we will store (so we don't have to store all)
print STDERR
    "\n\nreading clusterfile allwithdups.clusters to hash clustered ids (with added rotated repeats)...";
open my $cluster_fh, "<", "$rotatedfolder/allwithdups.clusters" or die $!;
$negcount = 0;
while (<$cluster_fh>) {

    #print STDERR $_;

    my @values = split( ',', $_ );

    foreach my $val (@values) {

        my $dir = '\'';
        if ( $val =~ /\"/ ) { $dir = '"'; }
        $val =~ s/[\'\"]//g;
        $val = trim($val);

        if ( $val > 0 ) {
            $RHASH{"$val"} = $dir;

            #print "setting `$val`\n";
        }
        else {
            $negcount++;
        }
    }

    #exit(1);
}
close($cluster_fh);

print STDERR
    keys(%RHASH)
    . " positive entries inserted into hash. (plus $negcount neg reference ones not in hash) ("
    . ( time() - $timestart )
    . ") secs\n\n";

# clear  tables
print STDERR "\ntruncating database tables\n\n";
my ( $trunc_query, $sth );
$dbh->begin_work;
$trunc_query = get_trunc_query( $run_conf{BACKEND}, "replnk" );
$sth         = $dbh->do($trunc_query);
$trunc_query = get_trunc_query( $run_conf{BACKEND}, "fasta_reads" );
$sth         = $dbh->do($trunc_query);
$dbh->commit;

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");
$dbh->do("PRAGMA journal_mode = TRUNCATE");

# Insert all ref TRs from replnk file
$sth = $dbh->prepare(
    qq{INSERT INTO
    replnk(rid,sid,first,last,patsize,copynum,pattern,profile,profilerc,profsize)
    VALUES (?,?,?,?,?,?,?,?,?,?)}
);

$timestart = time();
print STDERR
    "\nreading index file and storing relevant entries in database..."
    . "\n\n";
opendir( DIR, $indexfolder );
my @filelist   = readdir(DIR);
my @indexfiles = sort grep( /\.(?:index\.renumbered)$/, @filelist );
my @readfiles  = sort grep( /\.(?:reads)$/, @filelist );
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
my @replnk_rows;

die "read file count doesn't equal index file count: ($indexcount vs "
    . 1 * @readfiles . ")\n"
    if ( $indexcount != @readfiles );

foreach my $ifile (@indexfiles) {

    open( $fh1, "<", "$indexfolder/$ifile" ) or die $!;
    $i = 0;

    my $lfile = $ifile;
    $lfile =~ s/index/leb36/;

    open( $fh2, "<", "$indexfolder/$lfile" ) or die $!;

    print STDERR "\n" . $ifile . "-" . $lfile . "...";

    my $line1 = read_file_line($fh1);
    my $line2 = read_file_line($fh2);

    while ( $line1 && $line2 ) {
        if ( $line1
            =~ /^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d)\t(\d+)\t([A-Z]+)$/ )
        {

            if ( exists $RHASH{"$1"} ) {
                $processed++;

                $id      = $1;
                $head    = $2;
                $first   = $3;
                $last    = $4;
                $copy    = $5;
                $pat     = $6;
                $pattern = $7;
                $i++;

                $head = trim($head);

                unless ( exists $HEADHASH{"$head"} ) {

          #$sth0->execute("$head") or die "Cannot execute: " . $sth->errstr();
                    $HEADHASH{"$head"} = $processed;
                }

       #print $id." ".$head." ".$first." ".$last." ".$copy." ".$pat." "." \n";

                my @values = split( ' ', $line2 );
                if ( $values[0] != $id ) {
                    die
                        "id from index file ($id) does not match id from leb36 file ($values[0])";
                }

                my $profile   = $values[5];
                my $profilerc = $values[6];
                my $proflen   = length($profile) / 2;

#$sth->execute($id,$HEADHASH{"$head"},$first,$last,$copy,$pat,$pattern,$profile,$profilerc,$proflen) or die "Cannot execute: " . $sth->errstr();

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
                    if ($rows) {
                        $inserted += $rows;
                        @replnk_rows = ();
                    }
                    else {
                        die
                            "Something went wrong inserting, but somehow wasn't caught!\n";
                    }
                }

            }

            #if ($i>10) { last; }
        }

        #print $_;
        $line1 = read_file_line($fh1);
        $line2 = read_file_line($fh2);

    }
    close($fh1);
    close($fh2);
    print STDERR $i;

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
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

if ( $inserted == keys(%RHASH) ) {
    print STDERR "\n\n..."
        . $inserted
        . " read repeat entries inserted into database ("
        . ( time() - $timestart )
        . ") secs.\n\n";
}
else {
    die "\n\nERROR: hash contains " .
        keys(%RHASH)
        . " entries, while index files only have $processed matching entries and only $inserted were inserted into the database. Aborting!\n\n";
}

#goto AAA;

# get the read files
$totalReads = 0;
$inserted   = 0;
$processed  = 0;
$timestart  = time();
print STDERR
    "\nreading input read files and storing relevant entries in database..."
    . "\n\n";

$sth
    = $dbh->prepare(
    qq{INSERT INTO fasta_reads (sid, head, dna) VALUES(?, ?, ?)})
    or die "Couldn't prepare statement: " . $dbh->errstr;

my $headstr       = "";
my $dnastr        = "";
my $qualstr       = "";
my $HEADER_SUFFIX = "";
my @fasta_reads_rows;

# print Dumper(\%HEADHASH) . "\n";
my $files_processed = 0;
for my $read_file (@readfiles) {
    open my $r_fh, "<", "$indexfolder/$read_file";
    while ( my $line = <$r_fh> ) {
        chomp $line;
        ( $headstr, $dnastr ) = split "\t", $line;

        # Special last line
        if ( $headstr eq 'totalreads' ) {
            $totalReads += $dnastr;

            # Jump out of while loop
            last;
        }
        $headstr = trim($headstr);
        $dnastr  = trimall($dnastr);

        # warn "head: $headstr\tdna: $dnastr\n";

        my $dnabak = $dnastr;
        if ( exists $HEADHASH{"$headstr"} ) {

            $processed++;

            if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                warn "Read does not start with keyseq TCAG : "
                    . $dnabak . " ("
                    . $headstr . ")\n";
            }

            push @fasta_reads_rows,
                [ $HEADHASH{"$headstr"}, "$headstr", "$dnastr" ];

            if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                my $cb   = gen_exec_array_cb( \@fasta_reads_rows );
                my $rows = vs_db_insert( $dbh, $sth, $cb,
                    "Error inserting reads. row dump: "
                        . Dumper( \@fasta_reads_rows ) );
                if ($rows) {
                    $inserted += $rows;
                    @fasta_reads_rows = ();
                }
                else {
                    die
                        "Something went wrong inserting, but somehow wasn't caught!\n";
                }
            }
        }
    }

    close $r_fh;
    $files_processed++;
    say STDERR " (processed: $processed)";
}

# cleanup
if (@fasta_reads_rows) {
    my $cb   = gen_exec_array_cb( \@fasta_reads_rows );
    my $rows = vs_db_insert( $dbh, $sth, $cb,
        "Error inserting reads. row dump: " . Dumper( \@fasta_reads_rows ) );
    if ($rows) {
        $inserted += $rows;
        @fasta_reads_rows = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

# reenable indices
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");

# disconnect
$dbh->disconnect();

set_statistics( { NUMBER_READS => $totalReads } ) if ($totalReads);

# check
if ( $inserted == keys(%HEADHASH) ) {
    print STDERR "\n\n..."
        . $inserted
        . " reads inserted into database  ("
        . ( time() - $timestart )
        . ") secs..\n\n";
}
elsif ( $inserted != $processed ) {
    die "\n\nERROR: inserted into database "
        . $inserted
        . " entries, while input read files have $processed matching entries. Aborting!\n\n";
}
else {
    die "\n\nERROR: hash contains " .
        keys(%HEADHASH)
        . " entries, while input read files only have $inserted matching entries. Aborting!\n\n";
}

say STDERR "\n\nProcessing complete (insert_reads.pl).";

warn strftime( "\n\nend: %F %T\n\n", localtime );

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

# Perl function to remove whitespace from the start and end of the string
sub trim {
    my $string = shift;

    # Remove leading ">", if any
    $string =~ s/^>//;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Perl function to remove all whitespace
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
