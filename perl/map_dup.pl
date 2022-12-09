#!/usr/bin/env perl

# removes *reads* that are mapped to multiple references (due to multiple TRs). Picks the best TR by profile, when same, picks best by flank. When same, removes both.
# !!!  makes an exception if references are on same chromosome and close together

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh gen_exec_array_cb vs_db_insert set_statistics);

# Arguments
my $argc = @ARGV;
die "Usage: map_dup.pl expects 2 arguments.\n"
    unless $argc >= 2;

my $curdir   = getcwd();
my $cnf      = $ARGV[0];
my $outfile  = $ARGV[1];

my %run_conf = get_config("CONFIG", $cnf);
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

# TODO Better default or calculate in advance
my $maxRepeatsPerRead = 2;
my $RECORDS_PER_INFILE_INSERT = 100000;

$dbh->do("PRAGMA foreign_keys = OFF");

my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } );
my ($numTRsInRead) = $read_dbh->selectrow_array(q{
    SELECT count(*)
    FROM map JOIN replnk ON replnk.rid = map.readid
      JOIN rank USING (refid, readid)
      JOIN rankflank USING (refid, readid)
      JOIN fasta_ref_reps ON fasta_ref_reps.rid = map.refid
    WHERE replnk.sid = ? AND bbb = 1
    ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC})
    or die "Couldn't prepare statement: " . $dbh->errstr;

my $trsInRead_sth = $read_dbh->prepare(q{
    SELECT map.readid, map.refid, fasta_reads.head, rank.score, rankflank.score,
      reftab.head AS refhead, reftab.firstindex, reftab.lastindex, length(DNA)
    FROM map JOIN replnk ON replnk.rid = map.readid
      JOIN rank USING (refid, readid)
      JOIN rankflank USING (refid, readid)
      JOIN refdb.fasta_ref_reps reftab ON reftab.rid = map.refid
      JOIN fasta_reads using(sid)
    WHERE replnk.sid = ? AND bbb = 1
    ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC})
    or die "Couldn't prepare statement: " . $dbh->errstr;

# first get list of reads with multiple TRs that are mapped to more than one reference (sorted by read id)
my $deleted                            = 0;
my $ReadsDeleted                       = 0;
my $readsWithMultTRsMappedMultRefs_sth = $dbh->prepare(q{
    SELECT replnk.sid, count(map.refid)
    FROM map JOIN replnk ON replnk.rid = map.readid
    WHERE bbb = 1
    GROUP BY replnk.sid HAVING count(map.refid) > 1})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$readsWithMultTRsMappedMultRefs_sth->execute()
    or die "Cannot execute: " . $readsWithMultTRsMappedMultRefs_sth->errstr();


# create temp table for updates
$dbh->do(q{
CREATE TEMPORARY TABLE mduptemp (
    `refid` INT(11) NOT NULL,
    `readid` INT(11) NOT NULL,
    PRIMARY KEY (refid,readid) )})
    or die "Couldn't do statement: " . $dbh->errstr;
my $insert_mduptemp_sth = $dbh->prepare('INSERT INTO mduptemp VALUES(?, ?)');

my $i = 0;
my @mdup = ();
#open (my $outfh, ">$outfile");
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
        my $RefDiff = max( $oldRlast, $data2[7] ) - min( $oldRfirst, $data2[6] ) + 1;
        if ( $oldRfirst == -1 || $oldRlast == -1 ) { $RefDiff = 1000000;}

        #if ( $j == 1 ) {
        #    print $outfh ""
        #        . $i
        #        . ". read='"
        #        . $data2[2]
        #        . "' refs=("
        #        . $data[1] . ")\n";
        #}

# 1st entry can only be deleted due to NOT being on same chromosome and close together as 2nd entry (or more than $maxRepeatsPerRead entries exist)
# TODO Make this compare all adjacent pairs (not just TRs 1 and 2 in the read) for longer reads.
        if ($j == 2 && ( !( $data2[5] eq $oldrefhead && $RefDiff <= $readlen )
                         || ( $numTRsInRead > $maxRepeatsPerRead ))) {
            push @mdup, [ $oldrefid, $oldreadid ];
            #print $outfh "X";
            $deleted++;
            $isDeleted = 1;
        }

        #print $outfh "\n\t"
        #    . $data2[0] . "->"
        #    . $data2[1]
        #    . " P=$data2[3] F=$data2[4] ";

# delete every entry (except first which is deleted in another block) if more than $maxRepeatsPerRead exist
        if ( $j > 1 && $numTRsInRead > $maxRepeatsPerRead ) {
            push @mdup, [ $data2[1], $data2[0] ];
            #print $outfh "X";
            $deleted++;
            $isDeleted = 1;
        }

# else we are on 2nd etnry, delete 2nd entry if previous entry score is equal to this entry
#elsif ($j==2 && $data2[3]==$oldp && $data2[4]==$oldf) {
        elsif ( $j == $maxRepeatsPerRead ) {

# make an exception (DO NOT DELETE) if references are on same chromosome and close together as previous entry
            if ( $data2[5] eq $oldrefhead && $RefDiff <= $readlen ) { }
            else {
                push @mdup, [ $data2[1], $data2[0] ];
                #print $outfh "X";
                $deleted++;
                $isDeleted = 1;
            }
        }

        if ( @mdup % $RECORDS_PER_INFILE_INSERT == 0 ) {
            my $cb   = gen_exec_array_cb( \@mdup );
            my $rows = vs_db_insert( $dbh, $insert_mduptemp_sth, $cb,
                "Error when inserting entries into mduptemp table.\n" );
            @mdup = ();
        }

        $oldreadid  = $data2[0];
        $oldrefid   = $data2[1];
        $oldp       = $data2[3];
        $oldf       = $data2[4];
        $oldrefhead = $data2[5];
        $oldRfirst  = $data2[6];
        $oldRlast   = $data2[7];
    }
    $ReadsDeleted += $isDeleted;
}

# my $numReadsWithMultTRsMappedMultRefs
#     = $readsWithMultTRsMappedMultRefs_sth->rows;

if ( @mdup ) {
    my $cb   = gen_exec_array_cb( \@mdup );
    my $rows = vs_db_insert( $dbh, $insert_mduptemp_sth, $cb,
        "Error when inserting entries into mduptemp table.\n" );
    @mdup = ();
}


# update based on temp table
my $updfromtable = 0;
$updfromtable = $dbh->do(q{
    UPDATE map SET bbb=0
    WHERE EXISTS (
        SELECT refid FROM mduptemp t2
        WHERE map.refid = t2.refid AND map.readid = t2.readid)})
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->do("PRAGMA foreign_keys = ON");

if ( $updfromtable != $deleted ) {
    die "Updated bbb=0 number of entries ($updfromtable) not equal to"
        . " the number of deleted counter ($deleted), aborting!\n"
        . "You might need to rerun from step 12.";
}

# update BBB on stats
$dbh->do('UPDATE stats SET BBB = (SELECT count(*) FROM map WHERE bbb = 1)')
    or die "Couldn't do statement: " . $dbh->errstr;

#printf $outfh "%d entries deleted!\n", $deleted;
#printf $outfh "%d reads deleted!\n",   $ReadsDeleted;
#close $outfh;

$dbh->disconnect();
