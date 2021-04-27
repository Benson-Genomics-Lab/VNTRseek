#!/usr/bin/env perl

# removes *reads* that are mapped to multiple references (due to multiple TRs). Picks the best TR by profile, when same, picks best by flank. When same, removes both.
# !!!  makes an exception if references are on same chromosome and close together

my $RECORDS_PER_INFILE_INSERT = 100000;

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use POSIX qw(strftime);

use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh gen_exec_array_cb vs_db_insert set_statistics get_trunc_query);

print strftime( "\n\nstart: %F %T\n\n\n", localtime );

my $argc = @ARGV;

if ( $argc < 3 ) { die "Usage: map_dup.pl dbname msdir tempdir\n"; }

my $curdir = getcwd;

# TODO Better default or calculate in advance
my $maxRepeatsPerRead = 2;

my $DBSUFFIX = $ARGV[0];
my $run_dir  = $ARGV[1];
my $TEMPDIR  = $ARGV[2];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $DBSUFFIX, $run_dir );
my $dbh      = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

#my $sth3 = $dbh->prepare("UPDATE map SET bbb=0 WHERE refid=? AND readid=?;")
#                or die "Couldn't prepare statement: " . $dbh->errstr;

# create temp table for updates
my $query = q{CREATE TEMPORARY TABLE mduptemp (
    `refid` INT(11) NOT NULL,
    `readid` INT(11) NOT NULL,
    PRIMARY KEY (refid,readid)
    )};

$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

my $sth;
my $TEMPFILE;
$dbh->do("PRAGMA foreign_keys = OFF");

my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } );
my ($numTRsInRead) = $read_dbh->selectrow_array(
    q{SELECT COUNT(*)
    FROM map
    INNER JOIN replnk on replnk.rid=map.readid
    INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
    INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
    INNER JOIN fasta_ref_reps ON fasta_ref_reps.rid=map.refid
    WHERE replnk.sid=? AND bbb=1
    ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC}
) or die "Couldn't prepare statement: " . $dbh->errstr;

my $trsInRead_sth = $read_dbh->prepare(
    q{SELECT map.readid,map.refid,(
        SELECT head
        FROM fasta_reads
        WHERE fasta_reads.sid=replnk.sid
    ),rank.score,rankflank.score,reftab.head AS refhead,reftab.firstindex,reftab.lastindex,(
        SELECT length(DNA)
        FROM fasta_reads
        WHERE fasta_reads.sid=replnk.sid
    )
    FROM map
    INNER JOIN replnk on replnk.rid=map.readid
    INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
    INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
    INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid=map.refid
    WHERE replnk.sid=? AND bbb=1
    ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC}
) or die "Couldn't prepare statement: " . $dbh->errstr;

# first get list of reads with multiple TRs that are mapped to more than one reference (sorted by read id)
my $deleted                            = 0;
my $ReadsDeleted                       = 0;
my $readsWithMultTRsMappedMultRefs_sth = $dbh->prepare(
    q{SELECT replnk.sid,count(map.refid)
    FROM map
    INNER JOIN replnk on replnk.rid=map.readid
    WHERE bbb=1
    GROUP BY replnk.sid
    HAVING count(map.refid)>1}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$readsWithMultTRsMappedMultRefs_sth->execute()
    or die "Cannot execute: " . $readsWithMultTRsMappedMultRefs_sth->errstr();
my $insert_mduptemp_sth = $dbh->prepare(q{INSERT INTO mduptemp VALUES(?, ?)});

my $i = 0;
my @mdup = ();
while ( my @data = $readsWithMultTRsMappedMultRefs_sth->fetchrow_array() ) {

    # TODO Use read length to calculate $maxRepeatsPerRead
    $i++;

    $trsInRead_sth->execute( $data[0] );
    my $j          = 0;
    my $oldp       = -1;
    my $oldf       = -1;
    my $oldreadid  = -1;
    my $oldrefid   = -1;
    my $oldrefhead = "";
    my $oldRfirst  = -1;
    my $oldRlast   = -1;
    my $isDeleted  = 0;

    while ( my @data2 = $trsInRead_sth->fetchrow_array() ) {
        $j++;
        my $readlen = $data2[8];
        my $RefDiff
            = ( $oldRfirst == -1 || $oldRlast == -1 )
            ? 1000000
            : (
            max( $oldRlast, $data2[7] ) - min( $oldRfirst, $data2[6] ) + 1 );

        if ( $j == 1 ) {
            print "\n"
                . $i
                . ". read='"
                . $data2[2]
                . "' refs=("
                . $data[1] . ")";
        }

# 1st entry can only be deleted due to NOT being on same chromosome and close together as 2nd entry (or more than $maxRepeatsPerRead entries exist)
# TODO Make this compare all adjacent pairs (not just TRs 1 and 2 in the read) for longer reads.
        if ($j == 2
            && (  !( $data2[5] eq $oldrefhead && $RefDiff <= $readlen )
                || ( $numTRsInRead > $maxRepeatsPerRead ) )
            )
        {
            push @mdup, [ $oldrefid, $oldreadid ];
            print "X";
            $deleted++;
            $isDeleted = 1;
        }

        print "\n\t"
            . $data2[0] . "->"
            . $data2[1]
            . " P=$data2[3] F=$data2[4] ";

# delete every entry (except first which is deleted in another block) if more than $maxRepeatsPerRead exist
        if ( $j > 1 && $numTRsInRead > $maxRepeatsPerRead ) {

            push @mdup, [ $data2[1], $data2[0] ];
            print "X";
            $deleted++;
            $isDeleted = 1;
        }

# else we are on 2nd etnry, delete 2nd entry if previous entry score is equal to this entry
#elsif ($j==2 && $data2[3]==$oldp && $data2[4]==$oldf) {
        elsif ( $j == $maxRepeatsPerRead ) {

# make an exception (DO NOT DELETE) if references are on same chromosome and close together as previous entry
            if ( $data2[5] eq $oldrefhead && $RefDiff <= $readlen ) {

            }
            else {
                push @mdup, [ $data2[1], $data2[0] ];
                print "X";
                $deleted++;
                $isDeleted = 1;
            }
        }

        if ( ( @mdup % $RECORDS_PER_INFILE_INSERT == 0 ) ) {
            my $cb   = gen_exec_array_cb( \@mdup );
            my $rows = vs_db_insert( $dbh, $insert_mduptemp_sth, $cb,
                "Error when inserting entries into mduptemp table.\n" );
            if ($rows) {
                @mdup = ();
            }
            else {
                die
                    "Something went wrong inserting, but somehow wasn't caught!\n";
            }
        }

        $oldrefhead = $data2[5];
        $oldRfirst  = $data2[6];
        $oldRlast   = $data2[7];
        $oldp       = $data2[3];
        $oldf       = $data2[4];
        $oldreadid  = $data2[0];
        $oldrefid   = $data2[1];
    }

    $ReadsDeleted+=$isDeleted;

}

# my $numReadsWithMultTRsMappedMultRefs
#     = $readsWithMultTRsMappedMultRefs_sth->rows;

# $sth3->finish();
if ( @mdup ) {
    my $cb   = gen_exec_array_cb( \@mdup );
    my $rows = vs_db_insert( $dbh, $insert_mduptemp_sth, $cb,
        "Error when inserting entries into mduptemp table.\n" );
    if ($rows) {
        @mdup = ();
    }
    else {
        die
            "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}


# update based on temp table
my $updfromtable = 0;
$query = q{UPDATE map SET bbb=0
    WHERE EXISTS (
    SELECT refid FROM mduptemp t2
    WHERE map.refid = t2.refid AND map.readid=t2.readid
)};
$updfromtable = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->do("PRAGMA foreign_keys = ON");

if ( $updfromtable != $deleted ) {
    die
        "Updated bbb=0 number of entries($updfromtable) not equal to the number of deleted counter ($deleted), aborting! You might need to rerun from step 12.";
}

# update BBB on stats
$dbh->do(q{UPDATE stats SET BBB=(SELECT count(*) FROM map WHERE bbb=1)})
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->disconnect();

printf STDERR "\n\n%d entries deleted!\n", $deleted;
printf STDERR "%d reads deleted!\n",       $ReadsDeleted;

print STDERR strftime( "\n\nend: %F %T\n\n\n", localtime );

1;

