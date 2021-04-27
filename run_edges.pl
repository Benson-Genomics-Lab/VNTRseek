#!/usr/bin/env perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use List::Util qw[min max];

use strict;
use warnings;
use Cwd;
use DBI;
use POSIX qw(strftime);
use FindBin;
use File::Basename;
use File::Path qw(make_path);
use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics get_trunc_query);

sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $curdir = getcwd;

#my $PROCLU_PARAM = "$curdir/eucledian.dst $curdir/eucledian.dst 0 70 -se ";
my $PROCLU_PARAM = " $curdir/eucledian.dst 70 0 0 ";

my $files_to_process = 100;    # number of files to process in one batch
my $files_processed  = 0;      # files processed
my %p;                         # associates forked pids with output pipe pids

my $MYLOCK = 0;

my $argc = @ARGV;

if ( $argc < 8 ) {
    die
        "Usage: run_edges.pl reference_file edges_folder dbsuffix run_dir MINPROFSCORE NPROCESSORS PSEARCH TMPDIR\n";
}

my $inputfile    = $ARGV[0];
my $folder       = $ARGV[1];
my $DBSUFFIX     = $ARGV[2];
my $run_dir      = $ARGV[3];
my $MINPROFSCORE = $ARGV[4];
my $cpucount     = $ARGV[5];
my $PROCLU       = "$curdir/$ARGV[6]";
my $tmp          = $ARGV[7];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $DBSUFFIX, $run_dir );

my $clusters_processed = 0;
my $totalRefReps       = 0;
my $totalReadReps      = 0;

my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth2;
my $query;
my $result;
my $num;
my $i;
my $clusterid;
my $clusteridold;
my $repid;
my $exstring;
my $leftname  = "";
my $rightname = "";

# create folder
$exstring = "rm -f $folder -R";
system($exstring);
make_path($folder);

# Do away with tempmap
$dbh->do("CREATE TEMPORARY TABLE tempmap (rid int  PRIMARY KEY)")
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->do("INSERT INTO tempmap (rid) SELECT DISTINCT -refid FROM map")
    or die "Couldn't do statement: " . $dbh->errstr;

if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('ALTER TABLE tempmap DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('ALTER TABLE tempmap ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('SET AUTOCOMMIT = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('SET UNIQUE_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = OFF");

    # warn "\nTurning off AutoCommit\n";
    $dbh->{AutoCommit} = 0;
}

# print ref leb files
my %LHASH = ();
open( my $input_fh, "<$inputfile" . ".leb36" )
    or die "\nCannot open file '$inputfile'!\n";    # open for input
while (<$input_fh>) {                               # read file into list
    if (/^(\d+)/) {
        $LHASH{$1} = $_;
    }
}
close($input_fh);

$clusters_processed = 0;

# Do away with tempmap; use WHERE clusterlnk.repeatid IN
# (SELECT DISTINCT (-refid) FROM map)?
$query = qq{
    SELECT repeatid,clusterid
    FROM clusterlnk
        LEFT OUTER JOIN replnk ON clusterlnk.repeatid=replnk.rid
    WHERE clusterlnk.repeatid IN (SELECT DISTINCT (-refid) FROM map)
    ORDER BY clusterid
};
$sth = $dbh->prepare($query);
$sth->execute();
$i            = 0;
$clusteridold = -1;
open( my $outfile, ">$folder/dummy.txt" )
    or die "\nCannot open file '$folder/dummy.txt'!\n";    # open for output

while ( my @data = $sth->fetchrow_array() ) {
    $clusterid = $data[1];
    $repid     = $data[0];

    #print "repid: $repid"; exit(0);
    #print "$clusterid / $repid\n";

    if ( $clusterid != $clusteridold ) {
        $leftname = "refs.${clusterid}.leb36";
        close($outfile);
        open( $outfile, ">$folder/$leftname" )
            or die
            "\nCannot open file '$folder/$leftname'!\n";    # open for output

        $clusters_processed++;
        print "$clusters_processed\n";
    }

    {
        $repid = -1 * $repid;

        #$exstring = "echo -`grep $repid $inputfile`  >> $folder/$leftname";
        #system($exstring);

        if ( exists $LHASH{$repid} ) { print $outfile "-" . $LHASH{$repid}; }
        else {
            print STDERR "\nCould not lookup reference in hash ($repid)!\n";
            exit(1);
        }

    }

    $clusteridold = $clusterid;
    $i++;
}
close $outfile;

$num = $sth->rows;
$sth->finish;

print STDERR
    "Processing complete -- outputed $clusters_processed ref leb files. Creating input edges files ...\n\n";

# create edges files
$clusters_processed = 0;
$query              = qq{
    SELECT clusterid,refid,readid
    FROM map,clusterlnk
    WHERE -refid=clusterlnk.repeatid
    ORDER BY clusterid, refid, readid
    };
$sth = $dbh->prepare($query);
$sth->execute();
$i            = 0;
$clusteridold = -1;
my $refid;
my $readid;

while ( my @data = $sth->fetchrow_array() ) {
    $clusterid = $data[0];
    $refid     = $data[1];
    $readid    = $data[2];

    if ( $clusterid != $clusteridold ) {
        $leftname = "reads.${clusterid}.edgein";
        close(EFILE);
        open( EFILE, ">$folder/$leftname" )
            or die "Can't open $folder/$leftname!";
        $clusters_processed++;
        print "$clusters_processed\n";
    }
    print EFILE "-" . $refid . "," . $readid . "\n";

    $clusteridold = $clusterid;
    $i++;
}
$num = $sth->rows;
$sth->finish;

close(EFILE);

print STDERR
    "Processing complete -- outputed $clusters_processed edges input files.\n\n";

#AAA:

# $query = "DROP TABLE tempmap;";
# $dbh->do($query)
#     or die "Couldn't do statement: " . $dbh->errstr;

# $query = "CREATE TEMPORARY TABLE tempmap (rid int  PRIMARY KEY);";
# $dbh->do($query)
#     or die "Couldn't do statement: " . $dbh->errstr;

# if ( $run_conf{BACKEND} eq "mysql" ) {
#     $dbh->do('ALTER TABLE tempmap DISABLE KEYS;')
#         or die "Couldn't do statement: " . $dbh->errstr;
# }

# $query = "INSERT INTO tempmap(rid) SELECT DISTINCT readid FROM map;";
# $sth   = $dbh->do($query)
#     or die "Couldn't do statement: " . $dbh->errstr;

# if ( $run_conf{BACKEND} eq "mysql" ) {
#     $dbh->do('ALTER TABLE tempmap ENABLE KEYS;')
#         or die "Couldn't do statement: " . $dbh->errstr;
# }

# print read leb files
$clusters_processed = 0;

# $query              = "select cid FROM clusters ORDER BY cid;";
# $sth2               = $dbh->prepare($query);
# $sth2->execute();
$clusteridold = -1;

$query = qq{
    SELECT repeatid,clusterid,profile,profilerc,patsize,copynum
    FROM clusterlnk INNER JOIN replnk ON clusterlnk.repeatid=replnk.rid
        AND clusterlnk.repeatid IN (SELECT DISTINCT readid FROM map)
    ORDER BY clusterid
    };
$sth = $dbh->prepare($query);
$sth->execute();
my $myfile;
while ( my @data2 = $sth->fetchrow_array() ) {
    $clusterid = $data2[1];
    if ( $clusterid != $clusteridold ) {
        unless ( $clusteridold == -1 ) {
            close($myfile);
            $clusters_processed++;
            print "$clusters_processed\n";
        }

        # TODO Change this to not produce one file per cluster
        $rightname = "reads.${clusterid}.leb36";
        open( $myfile, ">$folder/$rightname" )
            or die "\nCannot open file '$folder/$rightname'!\n";
    }

#$query = "SELECT repeatid,clusterid,profile,profilerc,patsize,copynum FROM clusterlnk LEFT OUTER JOIN replnk ON clusterlnk.repeatid=replnk.rid INNER JOIN tempmap on clusterlnk.repeatid=tempmap.rid ORDER BY clusterid;";

    # $i = 0;
    # while ( my @data = $sth->fetchrow_array() ) {
    $repid = $data2[0];
    print "$clusterid / $repid\n";

    {
        my $copies = sprintf( "%.2lf", $data2[5] );
        my $pline
            = $repid . " "
            . $data2[4] . " "
            . $copies . " "
            . ( length( $data2[2] ) / 2 ) . " "
            . ( length( $data2[3] ) / 2 ) . " "
            . $data2[2] . " "
            . $data2[3]
            . " 0 0 0 0 |\n";
        print $myfile $pline;
    }

    # $i++;
    # }
    # $sth->finish;

    # close($myfile);
    # Shouldn't happen (?): if this cluster happens to not be in
    # the clusterlnk table, unlink the file we made.
    # if ( $sth->rows == 0 ) {
    #     unlink("$folder/$rightname");
    # }
    $clusteridold = $clusterid;
}
$sth->finish;
close($myfile);

print STDERR
    "Processing complete -- outputed $clusters_processed read leb files.\n\n";

$dbh->disconnect();

#AAA:

# get a list of input files
opendir( DIR, $folder );

# the only extensions are .leb36
my @tarballs = grep( /reads\.(\d+)\.(?:leb36)$/, readdir(DIR) );
closedir(DIR);
my $tarball_count = @tarballs;

if ( $tarball_count == 0 ) {
    warn strftime( "\n\nend: %F %T\n\n", localtime );
    warn "No files to process. Exiting...\n";
    exit 0;
}

print STDERR "$tarball_count supported files found in $folder\n";
$files_to_process = $tarball_count if $files_to_process > $tarball_count;

# enter dir
chdir($folder);

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $cpucount; $i++ ) { $p{ fork_proclu() } = 1 }

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {

    # check return value
    my ( $rc, $sig, $core ) = ( $? >> 8, $? & 127, $? & 128 );
    if ($core) {
        print STDERR "proclu process $pid dumped core\n";
        exit(1000);
    }
    elsif ( $sig == 9 ) {
        print STDERR "proclu process $pid was murdered!\n";
        exit(1001);
    }
    elsif ( $rc != 0 ) {
        print STDERR "proclu process $pid has returned $rc!\n";
        exit($rc);
    }

    if ( $p{$pid} ) {

        # one instance has finished processing -- start a new one
        delete $p{$pid};
        $p{ fork_proclu() } = 1;
    }
    else {
        die "************ Process $pid finished (not in hash)\n";
    }
}

print STDERR
    "Processing complete -- processed $files_processed cluster(s).\n";
warn strftime( "\n\nend: %F %T\n\n", localtime );

# update database
#AAA:

my %RHASH = ();
my %SHASH = ();

chdir($curdir);

$dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('SET AUTOCOMMIT = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('SET UNIQUE_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = OFF");

    # warn "\nTurning off AutoCommit\n";
    $dbh->{AutoCommit} = 0;
}

#$query = "UPDATE clusters SET profdensity=? WHERE cid=?;";
#$sth = $dbh->prepare($query);

# update cluster entries via temp file
print STDERR "Updating cluster table...\n";

my $TEMPFILE;
$query = q{
    CREATE TEMPORARY TABLE clustemp (
        `cid` integer NOT NULL PRIMARY KEY,
        `pd` float NOT NULL
    )};
if ( $run_conf{BACKEND} eq "mysql" ) {
    open( $TEMPFILE, ">$tmp/pcd_$DBSUFFIX.txt" ) or die $!;
    $query .= q{ENGINE=INNODB};
    $dbh->do($query)
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('ALTER TABLE clustemp DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth
        = $dbh->prepare(
        "LOAD DATA LOCAL INFILE '$tmp/pcd_$DBSUFFIX.txt' INTO TABLE clustemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do($query)
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->prepare(q{INSERT INTO clustemp VALUES(?, ?)})
        or die "Couldn't prepare statement: " . $dbh->errstr;
}

opendir( DIR, $folder );
@tarballs = grep( /(?:edges)$/, readdir(DIR) );
closedir(DIR);
$tarball_count = @tarballs;

$i = 0;
my $id     = "0";
my $ds     = 0.0;
my $pcd    = 0;
my $pcdupd = 0;
foreach (@tarballs) {
    $i++;

    #print $i.". ". $_ ."\n";
    if (/(\d+)\.leb36\.edges/i) {
        $id = $1;
    }

    open( my $fh, "$folder/$_" ) or die;
    while ( my $line = <$fh> ) {

        if ( $line =~ /-(\d+)(['"]),(\d+),(\d+\.\d+)/i ) {

            # do a check here if the read-ref pair is in the map file
            #$sth2->execute($1,$3);
            #$num = $sth2->rows;
            #if ($num>=1)
            {
                #print $1.$2." ".$3." ".$4."\n";
                $ds = $4;
                if ( $ds >= $MINPROFSCORE ) {

                    if ( exists $SHASH{$3} ) {
                        if ( $ds > $SHASH{$3} ) {
                            $RHASH{$3} = $1;
                            $SHASH{$3} = $ds;
                        }
                        elsif ( $ds == $SHASH{$3} ) {
                            $RHASH{$3} .= ( "," . $1 );
                            $SHASH{$3} = $ds;
                        }
                    }
                    else {
                        $RHASH{$3} = $1;
                        $SHASH{$3} = $ds;
                    }
                }
            }

        }
        elsif ( $line =~ /ave: (\d+\.\d+)/i ) {

            #print $id. ": ". $1."\n";
            #$sth->execute($1,$id);

            $pcd++;

            # warn "Inserting... $id, $i\n";
            if ( $run_conf{BACKEND} eq "mysql" ) {
                print $TEMPFILE "$id,$1\n";

                if ( $pcd % $RECORDS_PER_INFILE_INSERT == 0 ) {
                    close($TEMPFILE);
                    $pcdupd += $sth->execute();
                    open( $TEMPFILE, ">$tmp/pcd_$DBSUFFIX.txt" ) or die $!;
                }
            }
            elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                $sth->execute( $id, $1 );
                $pcdupd++;
            }
        }
    }
    close $fh;

}

# finish insert and update based on temp table
if ( $run_conf{BACKEND} eq "mysql" ) {
    close($TEMPFILE);
    $pcdupd += $sth->execute();
    $sth->finish;
    unlink("$tmp/pcd_$DBSUFFIX.txt");

    $sth = $dbh->do('ALTER TABLE clustemp ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}

$dbh->commit;

my $updfromtable = 0;

# $query = "UPDATE clustemp, clusters SET clusters.profdensity=clustemp.pd WHERE clusters.cid = clustemp.cid";
$query = q{
    UPDATE clusters SET profdensity=(
        SELECT pd FROM clustemp t2
        WHERE clusters.cid = t2.cid
    )
    WHERE EXISTS (
        SELECT * FROM clustemp t2
        WHERE clusters.cid = t2.cid
    )};
$updfromtable = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

# populate rank table
print STDERR "Populating rank table...\n";
$i = 0;
my $j       = 0;
my $rankins = 0;

# prepare rank table
$sth = $dbh->do( get_trunc_query( $run_conf{BACKEND}, "rank" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('ALTER TABLE rank DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
    open( $TEMPFILE, ">$tmp/ranktemp_$DBSUFFIX.txt" ) or die $!;

    $query
        = "LOAD DATA LOCAL INFILE '$tmp/ranktemp_$DBSUFFIX.txt' INTO TABLE rank FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';";
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $query = q{
        INSERT INTO rank VALUES(?, ?, ?, ?, "'")};
}

# insert rank entries via temp file
$sth = $dbh->prepare($query)
    or die "Couldn't prepare statement: " . $dbh->errstr;

while ( my ( $key, $value ) = each(%RHASH) ) {
    $i++;

    my @pieces = split( /,/, $value );
    my $ties = scalar(@pieces) - 1;
    if ( $ENV{DEBUG} ) {
        warn "$key => $value (" . $SHASH{$key} . "), ties: $ties\n";
    }
    foreach my $ps (@pieces) {

        $j++;

        if ( $run_conf{BACKEND} eq "mysql" ) {
            print $TEMPFILE $ps, ",", $key, ",", $SHASH{$key}, ",", $ties,
                "\n";
            if ( $j % $RECORDS_PER_INFILE_INSERT == 0 ) {
                close($TEMPFILE);
                $rankins += $sth->execute();
                open( $TEMPFILE, ">$tmp/ranktemp_$DBSUFFIX.txt" ) or die $!;
            }
        }
        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
            $sth->execute( $ps, $key, $SHASH{$key}, $ties );
            $rankins++;
        }
    }
}

# finish insert
if ( $run_conf{BACKEND} eq "mysql" ) {
    close($TEMPFILE);
    $rankins += $sth->execute();
}

$dbh->commit;
warn "Inserted ($j) rank records for $i reads.\n";

# create temp table for deletions
$query = q{CREATE TEMPORARY TABLE ranktemp (
    `refid` INTEGER NOT NULL,
    `readid` INTEGER NOT NULL,
    PRIMARY KEY (refid,readid)
)};

if ( $run_conf{BACKEND} eq "mysql" ) {
    $query .= ' ENGINE=INNODB';
}

$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

if ( $run_conf{BACKEND} eq "mysql" ) {

    # enable keys in rank, set stats, print msg
    $dbh->do('ALTER TABLE rank ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $dbh->do('ALTER TABLE ranktemp DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    open( $TEMPFILE, ">$tmp/ranktemp_$DBSUFFIX.txt" ) or die $!;
}

warn "Prunning (keep best ref TR for each read TR) from rank table...\n";
$query = q{Select refid,readid,sid,score
    FROM rank INNER JOIN
        replnk ON rank.readid=replnk.rid
    ORDER BY readid, score};
$sth = $dbh->prepare($query);
$sth->execute();
$i = 0;
my $count    = 0;
my $oldseq   = -1;
my $oldref   = -1;
my $oldread  = -1;
my $oldscore = -1.0;

if ( $run_conf{BACKEND} eq "sqlite" ) {
    $query = q{INSERT INTO ranktemp VALUES(?, ?)};
    $sth2  = $dbh->prepare($query);
}

while ( my @data = $sth->fetchrow_array() ) {
    if ( $data[1] == $oldread && $data[3] != $oldscore ) {

        # delete old one
        if ( $run_conf{BACKEND} eq "mysql" ) {
            print $TEMPFILE $oldref, ",", $oldread, "\n";
        }
        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
            $sth2->execute( $oldref, $oldread );
        }

        $count++;
    }
    $oldref   = $data[0];
    $oldread  = $data[1];
    $oldseq   = $data[2];
    $oldscore = $data[3];
    $i++;
}

$sth->finish;
$dbh->commit;
print STDERR "Prunning complete. Pruned $count rank records.\n";

print STDERR "Prunning (one TR/same read) from rank table...\n";

# readid added for tie resolution to keep rank and rankflank entries more in sync
$query = q{Select refid,readid,sid,score
    FROM rank INNER JOIN
        replnk ON rank.readid=replnk.rid
    ORDER BY refid, sid, score, readid};
$sth = $dbh->prepare($query);
$sth->execute();
$i       = 0;
$count   = 0;
$oldseq  = -1;
$oldref  = -1;
$oldread = -1;
while ( my @data = $sth->fetchrow_array() ) {

    if ( $data[0] == $oldref && $data[2] == $oldseq ) {

        # delete old one
        if ( $run_conf{BACKEND} eq "mysql" ) {
            print $TEMPFILE $oldref, ",", $oldread, "\n";
        }
        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
            $sth2->execute( $oldref, $oldread );
        }

        $count++;
    }
    $oldref  = $data[0];
    $oldread = $data[1];
    $oldseq  = $data[2];
    $i++;
}

$sth->finish;
$dbh->commit;
warn "Prunning complete. Pruned $count rank records.\n";

if ( $run_conf{BACKEND} eq "mysql" ) {

    # load the file into tempfile
    close($TEMPFILE);
    $sth
        = $dbh->do(
        "LOAD DATA LOCAL INFILE '$tmp/ranktemp_$DBSUFFIX.txt' INTO TABLE ranktemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
        ) or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('ALTER TABLE ranktemp ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    # delete from rank based on temptable entries
    $query = q{DELETE FROM t1
    USING rank t1 INNER JOIN
        ranktemp t2 ON (
            t1.refid = t2.refid AND t1.readid = t2.readid
        )};

    # cleanup temp file
    unlink("$tmp/ranktemp_$DBSUFFIX.txt");
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {

    # delete from rank based on temptable entries
    $query = qq{
        DELETE FROM rank
        WHERE EXISTS (
            SELECT * FROM ranktemp t2
            WHERE rank.refid = t2.refid
                AND rank.readid = t2.readid
        )};
}

my $delfromtable = 0;
$sth          = $dbh->prepare($query);
$delfromtable = $sth->execute();

# set old db settings
if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('SET AUTOCOMMIT = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET UNIQUE_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = ON");
    $dbh->{AutoCommit} = 1;
}

if ( $updfromtable != $pcdupd ) {
    die
        "Updated number of cluster entries($updfromtable) not equal to the number of inserted counter ($pcdupd), aborting! You might need to rerun from step 12.";
}
if ( $rankins != $j ) {
    die
        "Inserted number of rank entries($rankins) not equal to the number of inserted counter ($j), aborting! You might need to rerun from step 12.";
}
if ( $delfromtable != $count ) {
    die
        "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($count), aborting! You might need to rerun from step 12.";
}

$dbh->disconnect();
set_statistics(
    {   RANK_EDGES_OVERCUTOFF => $j,
        RANK_REMOVED_SAMEREF  => $count,
        RANK_REMOVED_SAMESEQ  => $count
    }
);

print STDERR "Finished. Deleted from rank using temptable: $delfromtable\n";
warn strftime( "\n\nend: %F %T\n\n", localtime );

1;

############################ Procedures ###############################################################

sub fork_proclu {
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
    print STDERR 'Processing files '
        . ( $files_processed + 1 ) . ' to '
        . ( $until + 1 ) . "\n";

    #my $output_prefix = "$root/$files_processed-$until";
    my @file_slice = @tarballs[ ($files_processed) .. ($until) ];
    my $file_slice_count = @file_slice;
    $files_processed += $file_slice_count;
    my $proclu_string;

    # unlock shared vars
    $MYLOCK = 0;

    defined( my $pid = fork )
        or die "Unable to fork: $!\n";

    # child
    if ( $pid == 0 ) {

        foreach (@file_slice) {

            #print STDERR "\t" . $_ . "\n";

            my $edgesin = $_;
            $edgesin =~ s/leb36/edgein/;

            my $reffile = $_;
            $reffile =~ s/read/ref/;

            #>/dev/null";
            #>$_.log";
            $proclu_string
                = $PROCLU . " "
                . $_ . " "
                . $reffile
                . $PROCLU_PARAM
                . "$edgesin > /dev/null";

#$proclu_string = $PROCLU . " " .  $_ . " "  . $reffile . $PROCLU_PARAM  . "$edgesin > ${edgesin}.proclu_log";
#print STDERR $proclu_string."\n";
#exit(1);
            system($proclu_string);
            if ( $? == -1 ) { die "command failed: $!\n"; }
            else {
                my $rc = ( $? >> 8 );

# empty clusters cause nonzero return code, cause some files missing
# if ( 0 != $rc ) { print "proclu returned $rc ( $proclu_string  )!"; exit($rc); }
            }

        }

        #print STDERR "\n";

        # child must never return
        exit 0;

        # parent
    }
    else {
        return $pid;
    }

    return 0;
}

