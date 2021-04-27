#!/usr/bin/env perl

# takes .index.seq files and calls pcr_dup program to remove duplicates

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use POSIX qw(strftime);
use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib";
use vutil
    qw(get_config get_dbh set_statistics get_trunc_query gen_exec_array_cb);

my $RECORDS_PER_INFILE_INSERT = 100000;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}
warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $argc = @ARGV;

if ( $argc < 7 ) {
    die
        "Usage: pcr_dup.pl indexfolder profcleanfolder dbname msdir cpucount tempdir ignorepcrdups\n";
}

my $curdir = getcwd;

my $indexfolder  = $ARGV[0];
my $pcleanfolder = $ARGV[1];
my $DBSUFFIX     = $ARGV[2];
my $run_dir      = $ARGV[3];
my $cpucount     = $ARGV[4];
my $TEMPDIR      = $ARGV[5];
my $KEEPPCRDUPS  = $ARGV[6];

# get run configuration, let's us connect to DB
my %run_conf = get_config( $DBSUFFIX, $run_dir );
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";
my %stats;

# process
warn "reading: $indexfolder\n";

opendir( my $dir, $indexfolder );
my @indexfiles = grep( /\.(?:seq)$/, readdir($dir) );
closedir($dir);

#my $i=0;
#foreach my $ifile (@indexfiles) {
# $i++;
# print "\n $i. ". $ifile."...";

# my $exstring = "./pcr_dup.exe $indexfolder/$ifile $indexfolder/$ifile.pcr_dup 2 3 > /dev/null";
# #my $exstring = "./pcr_dup.exe $indexfolder/$ifile $indexfolder/$ifile.pcr_dup 2 3 > $indexfolder/$ifile.pcr_log";
# system($exstring);
#}

my $files_to_process = 100;    # number of files to process in one batch
my $files_processed  = 0;      # files processed
my %p;                         # associates forked pids with output pipe pids

my $MYLOCK = 0;

my $tarball_count = @indexfiles;
warn "$tarball_count supported files found in $indexfolder\n";

#die "Exiting\n" if $tarball_count == 0;
$files_to_process = $tarball_count if $files_to_process > $tarball_count;

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $cpucount; $i++ ) { $p{ fork_pcrdup() } = 1 }

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {
    if ( $p{$pid} ) {

        # one instance has finished processing -- start a new one
        delete $p{$pid};
        $p{ fork_pcrdup() } = 1;
    }
    else {
        die "************ Process $pid finished (not in hash)\n";
    }
}

warn "Processing complete -- processed $files_processed cluster(s).\n";

# load results
warn "reading: $indexfolder\n";

# first count the intersect before pcr dup
my $rrintersect = 0;
my $sth         = $dbh->prepare(
    q{SELECT count(*)
    FROM rank INNER JOIN
        rankflank ON rank.refid=rankflank.refid AND rank.readid=rankflank.readid}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
($rrintersect) = $sth->fetchrow_array();
$sth->finish();
$stats{INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR} = $rrintersect;

# deleteing PCR DUPS in database
my %PENTRIES = ();

# create temp table for updates
my $query = q{CREATE TEMPORARY TABLE pduptemp (
    `refid` INT(11) NOT NULL,
    `readid` INT(11) NOT NULL,
    PRIMARY KEY (refid,readid)
    )};

$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");
$sth = $dbh->prepare(q{INSERT INTO pduptemp VALUES (?, ?)});

#my $query = "DELETE FROM rank WHERE refid=? AND readid=?;";
#$sth = $dbh->prepare($query);

#$query = "DELETE FROM rankflank WHERE refid=? AND readid=?;";
#my $sth1 = $dbh->prepare($query);

opendir( $dir, $indexfolder );
@indexfiles = grep( /\.(?:pcr_dup)$/, readdir($dir) );
closedir($dir);

my $i       = 0;
my $deleted = 0;
my @to_delete;
foreach my $ifile (@indexfiles) {

    $i++;

    if ( $ifile =~ /(\d+)\.seq\.pcr_dup/ ) {

        my $ref         = $1;
        my $filedeleted = 0;
        print "\n $i. " . $ifile . "...";

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

                    #my $read = max($1,$5);

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

                        #$sth->execute($ref,$read);
                        #$sth1->execute($ref,$read);
                        #print "deleting: -$ref -> $read\n";
                        $deleted++;
                        $filedeleted++;
                        push @to_delete, [ $ref, $read ];

                        if ( $deleted % $RECORDS_PER_INFILE_INSERT == 0 ) {
                            my $cb = gen_exec_array_cb( \@to_delete );
                            my $rows = vs_db_insert( $dbh, $sth, $cb,
                                "Error when inserting entries into temporary pcr duplicates table.\n");
                            if ($rows) {
                                @to_delete = ();
                            }
                            else {
                                die
                                    "Something went wrong inserting, but somehow wasn't caught!\n";
                            }
                        }

                  #print "$1 (newid:$NEWIDS{$1}) - $5 (newid: $NEWIDS{$5})\n";
                  #exit(1);

                        $PENTRIES{ $ref . "_" . $read } = 1;

                    }

                    $RHASH{$read} = 1;

                }
            }

            close($fh);
            print "($filedeleted deleted)";

        }

    }

}

#$sth->finish;
#$sth1->finish;

if (@to_delete) {
    my $cb = gen_exec_array_cb( \@to_delete );
    my $rows = vs_db_insert( $dbh, $sth, $cb,
        "Error when inserting entries into temporary pcr duplicates table.\n");
    if ($rows) {
        @to_delete = ();
    }
    else {
        die
            "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

# delete from rankdflank based on temptable entries
$query = q{DELETE FROM rank
WHERE EXISTS (
    SELECT * FROM pduptemp t2
    WHERE rank.refid = t2.refid
        AND rank.readid = t2.readid
)};

$dbh->begin_work();
my $delfromtable = $dbh->do($query);
$dbh->commit();

# cleanup temp file
if ( $run_conf{BACKEND} eq "mysql" ) {
    unlink("$TEMPDIR/pduptemp_$DBSUFFIX.txt");
}

print "\n\nProcessing complete (pcr_dup.pl), deleted $deleted duplicates.\n";
@stats{qw( RANK_REMOVED_PCRDUP RANKFLANK_REMOVED_PCRDUP )}
    = ( $deleted, $deleted );

# for accounting of pcr dups
warn "\nMaking a list of pcr_dup removed...\n";
if ( open( my $fh, ">$pcleanfolder/result/$DBSUFFIX.pcr_dup.txt" ) ) {
    $i = 0;
    for my $key ( sort keys %PENTRIES ) {
        $i++;
        print $fh "$i\t-" . $key . "\n";
    }
    warn "PCR_DUP list complete with $i removed entries.\n";
    close($fh);
}

# first count the intersect
$rrintersect = 0;
$sth         = $dbh->prepare(
    q{SELECT count(*) FROM rank
    INNER JOIN rankflank ON rank.refid=rankflank.refid AND rank.readid=rankflank.readid}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
($rrintersect) = $sth->fetchrow_array();
$sth->finish();
$stats{INTERSECT_RANK_AND_RANKFLANK} = $rrintersect;

# now exclude ties, mark in map table and record the number
warn "Updating BEST BEST BEST map entries...\n";

# clear all bbb entries
$dbh->do("UPDATE map SET bbb=0;")
    or die "Couldn't do statement: " . $dbh->errstr;

# clear all pduptemp entries
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "pduptemp" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

if ( $run_conf{BACKEND} eq "mysql" ) {
    $dbh->do('ALTER TABLE pduptemp DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}

# clear all pduptemp entries
$sth = $dbh->prepare(
    q{INSERT INTO pduptemp
    SELECT map.refid, map.readid FROM map
    INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
    INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
    WHERE rank.ties=0 OR rankflank.ties=0}
);
$dbh->begin_work();
$sth->execute();
$dbh->commit();

if ( $run_conf{BACKEND} eq "mysql" ) {
    $dbh->do('ALTER TABLE pduptemp ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}

$query = q{UPDATE map SET bbb=1
    WHERE EXISTS (
    SELECT refid FROM pduptemp t2
    WHERE map.refid = t2.refid AND map.readid=t2.readid
)};
$dbh->begin_work();
$i = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;
$dbh->commit();

$stats{BBB_WITH_MAP_DUPS} = $i;

# make a list of ties
warn "Making a list of ties (references)...\n";
if ( open( my $fh, ">$pcleanfolder/result/$DBSUFFIX.ties.txt" ) ) {
    my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } );
    $query = qq{SELECT map.refid, max(bbb) as mbb,
            (select head from refdb.fasta_ref_reps where rid=map.refid) as chr,
            (select firstindex from refdb.fasta_ref_reps where rid=map.refid) as tind
        FROM map INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
        INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
        GROUP BY map.refid HAVING mbb=0 ORDER BY chr, tind};
    $sth = $read_dbh->prepare($query);
    $sth->execute();
    $i = 0;
    while ( my @data = $sth->fetchrow_array() ) {
        $i++;
        print $fh "$i\t-"
            . $data[0] . "\t"
            . $data[2] . "\t"
            . $data[3] . "\n";
    }
    $sth->finish;
    close($fh);
}
warn "Ties list complete with $i removed references.\n";

# make a list of ties
warn "\nMaking a list of ties (entries)...\n";
if ( open( my $fh, ">$pcleanfolder/result/$DBSUFFIX.ties_entries.txt" ) ) {

    $query = q{SELECT map.refid, map.readid,rank.ties,rankflank.ties
    FROM map
    INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
    INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
    WHERE bbb=0 ORDER BY map.refid,map.readid};
    $sth = $dbh->prepare($query);
    $sth->execute();
    $i = 0;
    while ( my @data = $sth->fetchrow_array() ) {
        $i++;
        print $fh "$i\t-"
            . $data[0] . "\t"
            . $data[1] . "\t"
            . $data[2] . "\t"
            . $data[3] . "\n";
    }
    $sth->finish;
    close($fh);
}
warn "Ties list complete with $i removed entries.\n";

# set old settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");

if ( $delfromtable != $deleted ) {
    die
        "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($deleted), aborting! You might need to rerun from step 12.";
}

$dbh->disconnect();
warn strftime( "\n\nend: %F %T\n\n", localtime );
set_statistics( \%stats );

1;

############################ Procedures ###############################################################

sub fork_pcrdup {
    if ( $files_processed >= $tarball_count ) {
        return 0;
    }

    # wait for shared variables to unlock
    while ($MYLOCK) { }

    # lock shared vars
    $MYLOCK = 1;

    # use a predefined number of files
    my $until = $files_processed + $files_to_process - 1;
    $until = $tarball_count - 1 if $until > ( $tarball_count - 1 );
    warn 'Processing files '
        . ( $files_processed + 1 ) . ' to '
        . ( $until + 1 ) . "\n";

    #my $output_prefix = "$root/$files_processed-$until";
    my @file_slice = @indexfiles[ ($files_processed) .. ($until) ];
    my $file_slice_count = @file_slice;
    $files_processed += $file_slice_count;
    my $exstring;

    # unlock shared vars
    $MYLOCK = 0;

    defined( my $pid = fork )
        or die "Unable to fork: $!\n";

    # child
    if ( $pid == 0 ) {

        foreach (@file_slice) {

            #warn "\t" . $_ . "\n";

            my $exstring
                = "./pcr_dup.exe $indexfolder/${_} $indexfolder/${_}.pcr_dup 0 2 $KEEPPCRDUPS > /dev/null";
            system($exstring);

        }

        #warn "\n";

        # child must never return
        exit 0;

        # parent
    }
    else {
        return $pid;
    }

    return 0;
}

