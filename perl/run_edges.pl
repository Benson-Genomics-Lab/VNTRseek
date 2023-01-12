#!/usr/bin/env perl

use List::Util qw[min max];
use strict;
use warnings;
use Cwd;
use DBI;
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics);

# Arguments
die "Usage: run_edges.pl expects 6 arguments.\n"
    unless scalar @ARGV >= 6;

my $curdir       = getcwd();
my $inputfile    = $ARGV[0];
my $folder       = $ARGV[1];
my $cnf          = $ARGV[2];
my $MINPROFSCORE = $ARGV[3];
my $cpucount     = $ARGV[4];
my $PROCLU       = "$curdir/$ARGV[5]";


# load config, prep SQL driver, and set constants
my %run_conf = get_config("CONFIG", $cnf);
my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

my $PROCLU_PARAM = "$curdir/eucledian.dst 70 0 0";
my $files_to_process = 1000;

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");


# gather reference entries
my %LHASH = ();
open(my $input_fh, "<$inputfile" . ".leb36")
    or die "\nCannot open file '$inputfile'!\n";
while (<$input_fh>) {     # read file into list
    if (/^(\d+)/) {
        $LHASH{-$1} = $_; # note sign inversion, newlines not stripped!
    }
}
close($input_fh);


# create map for generating reference leb files.
my $query = q{
    SELECT clusterid, GROUP_CONCAT(repeatid, " ")
    FROM clusterlnk
        LEFT OUTER JOIN replnk ON clusterlnk.repeatid=replnk.rid
    WHERE clusterlnk.repeatid IN (SELECT DISTINCT (-refid) FROM map)
    GROUP BY clusterid;};
my $sth = $dbh->prepare($query);
$sth->execute();

my $clusterid;
my $repids;
my %refmap = ();
while (($clusterid, $repids) = $sth->fetchrow_array()) {
    $refmap{$clusterid} = $repids;
}
$sth->finish();
print "Reference map loaded.\n";


# create map for generating edges files
# x'0a' is a newline character
# sql does not convert \n
$query = q{
    SELECT clusterid, GROUP_CONCAT(datum, x'0a')
    FROM (SELECT clusterid, ("-" || refid || "," || readid) as datum
          FROM map, clusterlnk
          WHERE -refid=clusterlnk.repeatid
          ORDER BY refid, readid) AS t
    GROUP BY clusterid;};
$sth = $dbh->prepare($query);
$sth->execute();

my $text;
my %edgemap = ();
while (($clusterid, $text) = $sth->fetchrow_array()) {
    $edgemap{$clusterid} = $text;
}
$sth->finish();
print "Edge map loaded.\n";


# prep sample lebs
# x'0a' is a newline character
# sql does not convert \n
$query = q{
    SELECT clusterid, GROUP_CONCAT(datum, x'0a')
    FROM (SELECT clusterid, printf("%d %s %.2lf %d %d %s %s 0 0 0 0 |",
                                   repeatid, patsize, copynum,
                                   length(profile)/2, length(profilerc)/2,
                                   profile, profilerc) as datum
         FROM clusterlnk
             JOIN replnk ON clusterlnk.repeatid=replnk.rid
         WHERE clusterlnk.repeatid IN (SELECT DISTINCT readid FROM map)
         ORDER BY clusterid)
    GROUP BY clusterid};
$sth = $dbh->prepare($query);
$sth->execute();
print "Sample leb data ready.\n";


# move to working folder
chdir($folder);


# proclu forks
my $clusters_processed = 0;
my $coconut = -1;
my %p;
for (my $i = 0; $i < $cpucount; $i++) { $p{fork_proclu()} = 1;}

# wait for processes to finish and then fork new ones
while ((my $pid = wait) != -1) {

    # check return value
    my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
    if ($core) {
        warn "proclu process $pid dumped core\n";
        exit(1000);
    }
    elsif ($sig == 9) {
        warn "proclu process $pid was murdered!\n";
        exit(1001);
    }
    elsif ($rc != 0) {
        warn "proclu process $pid has returned $rc!\n";
        exit($rc);
    }

    # KA: not sure why we're checking this
    die "Unrelated process $pid finished?\n" unless $p{$pid};

    # one instance has finished processing -- start a new one
    delete $p{$pid};
    $p{fork_proclu()} = 1;
}
$sth->finish();
$coconut++; # since we used 0 indexing
print "Processing complete -- processed $clusters_processed cluster(s).\n";
chdir($curdir);


# update clusters tables
print "Updating clusters table.\n";

# make temp for clusters table
$query = q{
    CREATE TEMPORARY TABLE clustemp (
        `cid` integer NOT NULL PRIMARY KEY,
        `pd` float NOT NULL
    )};
$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;


# populate clustemp table and collect read scores while we're at it.
$sth = $dbh->prepare('INSERT INTO clustemp VALUES(?, ?)')
    or die "Couldn't prepare statement: " . $dbh->errstr;
my %RHASH = ();
my %SHASH = ();
my $pcdupd = 0;
for (my $c = 0; $c < $coconut; $c++) {

    # collect read scores
    # supposedly each read only gets assigned one cluster, and so only one score at this stage
    open(my $outfh, "$folder/scores.$c.txt") or die "Couldn't open file $folder/scores.$c.txt";
    while (my $line = <$outfh>) {
        my ($repid, $readids, $score) = split / /, $line;

        if (!exists $SHASH{$repid} || $score > $SHASH{$repid}) {
            $RHASH{$repid} = $readids;
            $SHASH{$repid} = $score;
        }
        elsif ($score == $SHASH{$3}) {
            $RHASH{$repid} .= ( "," . $readids );
        }
    }

    # collect cluster averages
    open($outfh, "$folder/aves.$c.txt") or die "Couldn't open file $folder/aves.$c.txt";
    while (my $line = <$outfh>) {
        my ($clusterid, $score) = split / /, $line;

        $sth->execute($clusterid, $score)
            or die "Error inserting $clusterid into clustemp table";
        $pcdupd++;
    }
    close $outfh;

    unlink "$folder/scores.$c.txt", "$folder/aves.$c.txt";
}
$sth->finish();


# do update
# UPDATE clusters SET profdensity=pd FROM clustemp WHERE clusters.cid = t2.cid
$query = q{
    UPDATE clusters SET profdensity=(
        SELECT pd FROM clustemp t2
        WHERE clusters.cid = t2.cid
    )
    WHERE EXISTS (
        SELECT * FROM clustemp t2
        WHERE clusters.cid = t2.cid
    )};
my $updfromtable = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;
print "Updated $updfromtable clusters from $pcdupd processed.\n";


# populate rank table
print "Populating rank table.\n";

# prepare rank table
$sth = $dbh->do("DELETE FROM rank")
    or die "Couldn't do statement: " . $dbh->errstr;

$query = q{INSERT INTO rank VALUES(?, ?, ?, ?, "'")};
$sth = $dbh->prepare($query)
    or die "Couldn't prepare statement: " . $dbh->errstr;

$dbh->begin_work();
my $nreads = 0;
my $rankins = 0;
while (my ($key, $value) = each(%RHASH)) {
    $nreads++;

    my @pieces = split /,/, $value;
    my $ties = scalar(@pieces) - 1;
    if ($ENV{DEBUG}) {
        warn "$key => $SHASH{$key}; $value, ties: $ties\n";
    }
    for my $ps (@pieces) {
        $sth->execute($ps, $key, $SHASH{$key}, $ties)
            or die "Error inserting ($ps, $key) into rank table";
        $rankins++;
    }
}
$sth->finish();
$dbh->commit();
print "Inserted $rankins rank records for $nreads reads.\n";


# create temp table for deletions
$query = q{
    CREATE TEMPORARY TABLE ranktemp (
        `refid` INTEGER NOT NULL,
        `readid` INTEGER NOT NULL,
        PRIMARY KEY (refid, readid)
    )};
$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;
my $sth2  = $dbh->prepare(q{INSERT INTO ranktemp VALUES(?, ?)});


# prune reference ids with worse scores
print "Prunning (keep best ref TR for each read TR) from rank table.\n";
$query = q{
    SELECT refid, readid, sid, score
    FROM rank INNER JOIN
        replnk ON rank.readid=replnk.rid
    ORDER BY readid, score};
$sth = $dbh->prepare($query);
$sth->execute();

my $count    = 0;
my $oldref   = -1;
my $oldread  = -1;
my $oldscore = -1.0;
while (my @data = $sth->fetchrow_array()) {
    if ($data[1] == $oldread && $data[3] != $oldscore) {
        # record old one (smaller) for deletion
        $sth2->execute($oldref, $oldread)
            or die "Error inserting ($oldref, $oldread) into ranktemp table";
        $count++;
    }
    $oldref   = $data[0];
    $oldread  = $data[1];
    $oldscore = $data[3];
}
$sth->finish();
print "Prunning complete. Pruned $count rank records.\n";


# prune duplicates?
print "Prunning (one TR/same read) from rank table.\n";
$query = q{
    SELECT refid, readid, sid, score
    FROM rank INNER JOIN
        replnk ON rank.readid=replnk.rid
    ORDER BY refid, sid, score, readid};
# readid added for tie resolution to keep rank and rankflank entries more in sync
$sth = $dbh->prepare($query);
$sth->execute();

$count   = 0;
$oldref  = -1;
$oldread = -1;
my $oldseq  = -1;
while (my @data = $sth->fetchrow_array()) {
    if ($data[0] == $oldref && $data[2] == $oldseq) {
        # record old one (identical?) for deletion
        $sth2->execute($oldref, $oldread);
        $count++;
    }
    $oldref  = $data[0];
    $oldread = $data[1];
    $oldseq  = $data[2];
}
$sth->finish();
print "Prunning complete. Pruned $count rank records.\n";


# delete from rank based on ranktemp entries
$query = qq{
    DELETE FROM rank
    WHERE EXISTS (
        SELECT * FROM ranktemp t2
        WHERE rank.refid = t2.refid
            AND rank.readid = t2.readid
    )};
my $delfromtable = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;
print "Deleted from rank using temptable: $delfromtable.\n";


# set old db settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");
$dbh->disconnect();


if ( $updfromtable != $pcdupd ) {
    die "Updated number of cluster entries($updfromtable) not equal to" .
        " the number of inserted counter ($pcdupd), aborting!".
        " You might need to rerun from step 12.";
}
if ( $delfromtable != $count ) {
    die "Deleted number of entries($delfromtable) not equal to" .
        " the number of deleted counter ($count), aborting!" .
        " You might need to rerun from step 12.";
}

set_statistics({
    RANK_EDGES_OVERCUTOFF => $rankins,
    RANK_REMOVED_SAMEREF  => $count,
    RANK_REMOVED_SAMESEQ  => $count
});

1;

############################ Procedures ###############################################################

sub fork_proclu {

    my @clusterids = ();
    my @row;
    my $outfh;
    my $datum;
    for (my $k = 0; $k < $files_to_process; $k++) {
        @row = $sth->fetchrow_array(); # one row to one clusterid
        if (@row == 0) { last;} # empty list means no more rows, means no more clusters

        ($clusterid, $datum) = @row;
        push @clusterids, $clusterid; # spool cluster ids for child to process

        my $readfile = "reads.$clusterid.leb36"; # output read data from sql handle to file
        open $outfh, ">", $readfile or die "Cannot open file '$readfile'!\n";
        print $outfh $datum;
        close($outfh);
        $clusters_processed++;
    }

    if (scalar @clusterids == 0) { return 0;} # none spooled means no more clusters in feed

    $coconut++; # only increment if we know we're going to fork, vis it's the next line
    defined(my $pid = fork) or die "Unable to fork: $!\n";
    if ($pid != 0) { return $pid;} # ← parent, child ↓

    open(my $avesfh, ">aves.$coconut.txt")
      or die "Couldn't open output aves.$coconut.text";
    open(my $scoresfh, ">scores.$coconut.txt")
      or die "Couldn't open output scores.$coconut.text";

    for $clusterid (@clusterids) {
        # file names for proclu
        my $readfile = "reads.$clusterid.leb36";
        my $reffile  = "refs.$clusterid.leb36";
        my $edgefile = "reads.$clusterid.edgein";

        # generate reference file
        open ($outfh, ">", $reffile)
          or die "Cannot open file '$reffile'!\n";
        for my $repid ((split / /, $refmap{$clusterid})) {
            #if (!exists $LHASH{$repid}) { die "Couldn't find repid $repid in reference map (LHASH).";}
            print $outfh "-" . $LHASH{$repid}; # newlines not stripped during LHASH loading
        }

        # generate edges file
        open ($outfh, ">", $edgefile)
          or die "Cannot open file '$edgefile'!\n";
        #if (!exists $edgemap{$clusterid}) { die "Couldn't find clusterid $clusterid in edge map.";}
        print $outfh $edgemap{$clusterid};
        close($outfh);

        # call proclu
        my $proclu_string = "$PROCLU $readfile $reffile $PROCLU_PARAM $edgefile > /dev/null";
        system($proclu_string);
        if ($? == -1) { die "Proclu failed on cluster $clusterid: $!\n";}
        else {
            my $rc = ($? >> 8);
            if ( 0 != $rc ) {
                die "Proclu returned error code $rc: $?\n";
            }
        } # I don't think this is used

        # add output to intermediate files.
        open($outfh, "<", "$reffile.edges")
          or die "Couldn't find proclu results at '$reffile.edges'.";

        my %rHASH;   # rHASH: ref.repid → read.repid(s)
        my %sHASH;   # sHASH: ref.repid → sup {reads} profscore
        while (my $line = <$outfh>) {
            if ($line =~ /-(\d+)['"],(\d+),(\d+\.\d+)/i and $3 >= $MINPROFSCORE) {
                if (!exists $sHASH{$2} or $3 > $sHASH{$2}) {
                    $rHASH{$2} = $1;
                    $sHASH{$2} = $3;
                }
                elsif ($3 == $sHASH{$2}) { $rHASH{$2} .= "," . $1;}
            }
            elsif ($line =~ /ave: (\d+\.\d+)/i) {
                print $avesfh $clusterid . " " . $1 . "\n";
            }
        }

        for my $repid (keys %rHASH) {
            print $scoresfh "$repid $rHASH{$repid} $sHASH{$repid}\n";
        }

        unlink $readfile, $reffile, $edgefile, "$reffile.edges";
    }
    exit 0; # child does not return
}
