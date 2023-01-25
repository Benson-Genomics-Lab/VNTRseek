#!/usr/bin/env perl

# takes .index.seq files and calls pcr_dup program to remove duplicates

use strict;
use warnings;
use Cwd;
use DBI;
use File::Basename;
use List::Util qw(min max);

use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics gen_exec_array_cb);

my $argc = @ARGV;
die "Usage: pcr_dup.pl expects 6 arguments.\n"
    unless $argc >= 6;

# Arguments
my $curdir       = getcwd();
my $indexfolder  = $ARGV[0];
my $out_folder   = $ARGV[1];
my $RUN_NAME     = $ARGV[2];
my $cnf          = $ARGV[3];
my $cpucount     = $ARGV[4];
my $KEEPPCRDUPS  = $ARGV[5];

my %run_conf = get_config("CONFIG", $cnf);
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";
my %stats;

my $RECORDS_PER_INFILE_INSERT = 100000;
my $files_to_process = 100;    # number of files to process in one batch


### Indexing
my $sql_clause = q{
  FROM map JOIN rank using(refid, readid)
      JOIN rankflank using(refid, readid)
      JOIN replnk ON replnk.rid=map.readid
      JOIN fasta_reads on fasta_reads.sid=replnk.sid
  ORDER BY map.refid, map.readid};
my ($num) = $dbh->selectrow_array( q{SELECT COUNT(*) } . $sql_clause )
    or die "Couldn't execute statement: " . $dbh->errstr;
my $sth = $dbh->prepare(
    q{SELECT map.refid, map.readid, replnk.sid, replnk.first, replnk.last,
      replnk.copynum, replnk.patsize, replnk.pattern, fasta_reads.dna}
        . $sql_clause )
    or die "Couldn't prepare statement: " . $dbh->errstr;

$sth->execute() or die "Couldn't execute: " . $sth->errstr;

print "Best best best records: $num\n";

my $oldref  = -1;
my $i = 0;
my $nrefs = 0;

my $fn;
my $fh;
my @refindices = ();
while ( my @data = $sth->fetchrow_array() ) {
    if ( $data[0] != $oldref ) {
        if ( $i != 0 ) { close($fh); }
        $nrefs++; # Counts the number of .seq files
        open $fh, ">$indexfolder/$data[0].seq" or die $!;
        push @refindices, $data[0];
    }
    print $fh "$data[1] $data[2] $data[3] $data[4]";
    printf $fh " %.2lf", $data[5];
    $data[8] =~ s/\s+//g;
    print $fh " $data[6] $data[7] $data[8]\n";
    $oldref  = $data[0];
    $i++;
}
close ($fh) if ($fh);
$sth->finish();
$dbh->disconnect();

print "Indexing complete, $nrefs files created.\n";


### Calculating and Removing
my $numfiles = @refindices;
$files_to_process = $numfiles if $files_to_process > $numfiles;

my $files_processed  = 0;      # files processed
my %p;                         # forked pids

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $cpucount; $i++ ) { $p{ fork_pcrdup() } = 1 }

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {
    # check return value
    my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
    if ($core) {
        warn "pcr_dup process $pid dumped core\n";
        exit(1000);
    }
    elsif ($sig == 9) {
        warn "pcr_dup process $pid was murdered!\n";
        exit(1001);
    }
    elsif ($rc != 0) {
        warn "pcr_dup process $pid has returned $rc!\n";
        exit($rc);
    }

    die "Unrelated process $pid finished?\n" unless $p{$pid};

    # one instance has finished processing -- start a new one
    delete $p{$pid};
    $p{ fork_pcrdup() } = 1;
}
print "Processing complete -- processed $files_processed cluster(s).\n";


# first count the intersect before pcr dup
$sth = $dbh->prepare(q{SELECT count(*)
    FROM rank JOIN rankflank using(refid, readid)})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
($stats{INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR}) = $sth->fetchrow_array();
$sth->finish();

# deleteing PCR DUPS in database
my %PENTRIES = ();

# create temp table for updates
my $query = q{
  CREATE TEMPORARY TABLE pduptemp (
    `refid` INT(11) NOT NULL,
    `readid` INT(11) NOT NULL,
    PRIMARY KEY (refid,readid))};
$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;
$sth = $dbh->prepare(q{INSERT INTO pduptemp VALUES (?, ?)});

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");


# load results
print "Processing pcr_dup.exe output\n";

$i          = 0;
my $deleted = 0;
my @to_delete;
foreach my $ref (@refindices) {
    my $ifile = "$ref.seq.pcr_dup";

    $i++;
    my $filedeleted = 0;

    if ( open( my $fh, "$indexfolder/$ifile" ) ) {

        my %RHASH = ();

        # added at 1.02, to eliminate most connected nodes preferentiably
        my %RCOUNTS = ();
        my %NEWIDS  = ();
        while (<$fh>) {
            if (/^compare: (\d+) (\d+) (\d+) (\d+)\|(\d+)/) {
                $RCOUNTS{$1}++;
                $RCOUNTS{$5}++;
            }
        }
        my @keys = sort { $RCOUNTS{$a} <=> $RCOUNTS{$b} or $a <=> $b }
            keys %RCOUNTS;
        my $iditer = 1;
        foreach my $key (@keys) { $NEWIDS{$key} = $iditer; $iditer++; }

        seek $fh, 0, 0;
        while (<$fh>) {
            if (/^compare: (\d+) (\d+) (\d+) (\d+)\|(\d+)/) {

                my $read;
                if ( $NEWIDS{$1} > $NEWIDS{$5} ) {
                    $read = $1;
                }
                elsif ( $NEWIDS{$1} < $NEWIDS{$5} ) {
                    $read = $5;
                }
                else {
                    $read = max( $1, $5 );
                }

                if ( !exists $RHASH{$read} ) {
                    $deleted++;
                    $filedeleted++;
                    push @to_delete, [ $ref, $read ];

                    if ( $deleted % $RECORDS_PER_INFILE_INSERT == 0 ) {
                        my $cb = gen_exec_array_cb( \@to_delete );
                        my $rows = vs_db_insert( $dbh, $sth, $cb,
                            "Error when inserting entries into temporary pcr duplicates table.\n");
                        @to_delete = ();
                    }
                    $PENTRIES{ $ref . "_" . $read } = 1;
                }
                $RHASH{$read} = 1;
            }
        }
        close($fh);
    }
}

if (@to_delete) {
    my $cb = gen_exec_array_cb( \@to_delete );
    my $rows = vs_db_insert( $dbh, $sth, $cb,
        "Error when inserting entries into temporary pcr duplicates table.\n");
    @to_delete = ();
}
$sth->finish();

# delete from rank based on temptable entries
$query = q{
    DELETE FROM rank
    WHERE EXISTS (
        SELECT * FROM pduptemp t2
        WHERE rank.refid = t2.refid
          AND rank.readid = t2.readid)};
my $delfromtable = $dbh->do($query);


print "Processing complete, deleted $deleted duplicates.\n";
@stats{qw( RANK_REMOVED_PCRDUP RANKFLANK_REMOVED_PCRDUP )}
    = ( $deleted, $deleted );

# for accounting of pcr dups
# KA: turned off
#print "Making a list of pcr_dup removed.\n";
#if ( open( my $fh, ">$out_folder/$RUN_NAME.pcr_dup.txt" ) ) {
#    $i = 0;
#    for my $key ( sort keys %PENTRIES ) {
#        $i++;
#        print $fh "$i\t-" . $key . "\n";
#    }
#    print "PCR_DUP list complete with $i removed entries.\n";
#    close($fh);
#}

# re-count the intersect
$sth = $dbh->prepare(q{
    SELECT count(*)
    FROM rank JOIN rankflank using (refid, readid)})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
($stats{INTERSECT_RANK_AND_RANKFLANK}) = $sth->fetchrow_array();
$sth->finish();

# now exclude ties, mark in map table and record the number
print "Updating BEST BEST BEST map entries.\n";

# clear all bbb entries
$dbh->do("UPDATE map SET bbb=0;")
    or die "Couldn't do statement: " . $dbh->errstr;

# clear all pduptemp entries
$dbh->do( "DELETE FROM pduptemp" )
    or die "Couldn't do statement: " . $dbh->errstr;

# repopulate pduptemp
$sth = $dbh->prepare(q{
    INSERT INTO pduptemp
    SELECT map.refid, map.readid FROM map
      JOIN rank using(refid, readid)
      JOIN rankflank using(refid, readid)
    WHERE rank.ties = 0 OR rankflank.ties = 0});
$sth->execute();

$query = q{
    UPDATE map SET bbb=1
    WHERE EXISTS (
        SELECT refid FROM pduptemp
        WHERE map.refid = pduptemp.refid
          AND map.readid = pduptemp.readid)};
$stats{BBB_WITH_MAP_DUPS} = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;
$dbh->commit();

# make a list of ties
# KA: turned off
#print "Making a list of ties (references).\n";
#if ( open( my $fh, ">$out_folder/$RUN_NAME.ties.txt" ) ) {
#    my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } );
#    $query = q{
#        SELECT map.refid, max(bbb) as mbb,
#            (select head from refdb.fasta_ref_reps where rid=map.refid) as chr,
#            (select firstindex from refdb.fasta_ref_reps where rid=map.refid) as tind
#        FROM map JOIN rank using(refid, readid)
#          JOIN rankflank using(refid, readid)
#        GROUP BY map.refid HAVING mbb=0 ORDER BY chr, tind};
#    $sth = $read_dbh->prepare($query);
#    $sth->execute();
#    $i = 0;
#    while ( my @data = $sth->fetchrow_array() ) {
#        $i++;
#        print $fh "$i\t-"
#            . $data[0] . "\t"
#            . $data[2] . "\t"
#            . $data[3] . "\n";
#    }
#    $sth->finish();
#    close($fh);
#}
#print "Ties list complete with $i removed references.\n";

# make a list of ties
# KA: turned off
#print "Making a list of ties (entries).\n";
#if ( open( my $fh, ">$out_folder/$RUN_NAME.ties_entries.txt" ) ) {
#    $query = q{SELECT map.refid, map.readid, rank.ties, rankflank.ties
#    FROM map JOIN rank using(refid, readid)
#      JOIN rankflank using(refid, readid)
#    WHERE bbb=0 ORDER BY map.refid, map.readid};
#    $sth = $dbh->prepare($query);
#    $sth->execute();
#    $i = 0;
#    while ( my @data = $sth->fetchrow_array() ) {
#        $i++;
#        print $fh "$i\t-"
#            . $data[0] . "\t"
#            . $data[1] . "\t"
#            . $data[2] . "\t"
#            . $data[3] . "\n";
#    }
#    $sth->finish();
#    close($fh);
#}
#print "Ties list complete with $i removed entries.\n";

# set old settings
$dbh->commit();
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");
$dbh->disconnect();

if ( $delfromtable != $deleted ) {
    die "Deleted number of entries($delfromtable) not equal to the number of deleted"
        . " counter ($deleted), aborting! You might need to rerun from step 12.";
}

set_statistics( \%stats );

1;

############################ Procedures ###############################################################

sub fork_pcrdup {
    if ( $files_processed >= $numfiles ) { return 0;}

    # use a predefined number of files
    my $until = $files_processed + $files_to_process - 1;
    $until = $numfiles - 1 if $until > ( $numfiles - 1 );

    #my $output_prefix = "$root/$files_processed-$until";
    my @file_slice = @refindices[ ($files_processed) .. ($until) ];
    my $file_slice_count = @file_slice;
    $files_processed += $file_slice_count;

    defined( my $pid = fork )
        or die "Unable to fork: $!\n";

    if ( $pid != 0 ) { return $pid;} # parent
    # child
    foreach (@file_slice) {
        system("./pcr_dup.exe $indexfolder/${_}.seq $indexfolder/${_}.seq.pcr_dup"
            . " 0 2 $KEEPPCRDUPS > /dev/null");
    }
    # child must never return
    exit 0;
}
