#!/usr/bin/env perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use strict;
use warnings;
use List::Util qw[min max];
use Cwd;
use POSIX qw(strftime);
use DBI;

use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib";

use vutil
    qw(get_config get_dbh set_statistics get_trunc_query gen_exec_array_cb vs_db_insert);

sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $curdir = getcwd;

my $argc = @ARGV;
if ( $argc < 4 ) {
    die "Usage: run_flankcomp.pl inputfile dbname run_dir tempdir\n";
}

my $inputfile = $ARGV[0];
my $DBSUFFIX  = $ARGV[1];
my $run_dir   = $ARGV[2];
my $TEMPDIR   = $ARGV[3];

# get run config
my %run_conf = get_config( $DBSUFFIX, $run_dir );

my $clusters_processed = 0;
my $totalRefReps       = 0;
my $totalReadReps      = 0;

my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } )
    or die "Could not connect to database for reading: $DBI::errstr";
my $write_dbh = get_dbh()
    or die "Could not connect to database for writing: $DBI::errstr";

my $sth;
my $sth1;
my $sth2;
my $sth3;

my $mostReps    = 0;
my $mostRefReps = 0;
my $maxRange    = 0;

my $BREAK_SIZE = 4000;

#for my $i (keys %BELONG) {
#  print STDERR $i.":".$BELONG{$i}."\n";
#}
#exit(1);

# clear database cluster tables
$write_dbh->do( get_trunc_query( $run_conf{BACKEND}, "clusters" ) )
    or die "Couldn't do statement: " . $write_dbh->errstr;

$write_dbh->do( get_trunc_query( $run_conf{BACKEND}, "clusterlnk" ) )
    or die "Couldn't do statement: " . $write_dbh->errstr;

$write_dbh->do("PRAGMA foreign_keys = OFF");
$write_dbh->do("PRAGMA synchronous = OFF");

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

print STDERR "\nInserting into clusterlnk table...\n";

open my $fh, "<$inputfile" or die $!;
my $TEMPFILE;

$sth = $write_dbh->prepare(q{INSERT INTO clusterlnk VALUES (?, ?, ?, 0, 0)})
    or die "Couldn't prepare statement: " . $write_dbh->errstr;

# insert into clusterlnk
my $totalreps = 0;

# DBI->trace("3|SQL", "dbitrace.log");
my @cluster_links;
while (<$fh>) {
    $clusters_processed++;

    chomp;

    my @values = split( ',', $_ );

    foreach my $val (@values) {

        # insert clusterlnk entry

        $totalreps++;

        my $dir = q{'};
        if ( $val =~ s/([\'\"])// ) { $dir = $1; }
        push @cluster_links, [ $clusters_processed, $val, $dir ];

        if ( ( $totalreps % $RECORDS_PER_INFILE_INSERT == 0 ) ) {
            my $cb   = gen_exec_array_cb( \@cluster_links );
            my $rows = vs_db_insert( $write_dbh, $sth, $cb,
                "Error when inserting entries into our clusterlnk table.\n" );
            if ($rows) {
                @cluster_links = ();
            }
            else {
                die
                    "Something went wrong inserting, but somehow wasn't caught!\n";
            }
        }

    }

}    # end of while loop

# Finish insert
if (@cluster_links) {
    my $cb   = gen_exec_array_cb( \@cluster_links );
    my $rows = vs_db_insert( $write_dbh, $sth, $cb,
        "Error when inserting entries into our clusterlnk table.\n" );
    if ($rows) {
        @cluster_links = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

# Sync the DB so the next SELECTs work
$write_dbh->do("PRAGMA synchronous = FULL");
$write_dbh->do("PRAGMA synchronous = OFF");

print STDERR "\nPrinting DNA and inserting into cluster table...\n";

# now print dna and quals (also insert into cluster table)
$sth = $read_dbh->prepare(
    q{SELECT rid,flankleft,sequence,flankright,pattern,copynum,direction
    FROM refdb.fasta_ref_reps
    INNER JOIN clusterlnk ON rid=-repeatid
    WHERE clusterid = ?}
) or die "Couldn't prepare statement: " . $read_dbh->errstr;
$sth1 = $read_dbh->prepare(
    q{SELECT rid, dna, first, last, pattern, copynum,direction
    FROM fasta_reads
    INNER JOIN replnk ON fasta_reads.sid=replnk.sid
    INNER JOIN clusterlnk ON rid=repeatid
    WHERE clusterid = ?}
) or die "Couldn't prepare statement: " . $read_dbh->errstr;
$sth2 = $write_dbh->prepare(
    q{INSERT INTO clusters(cid,minpat,maxpat,repeatcount,refcount)
    VALUES(?,?,?,?,?)}
) or die "Couldn't prepare statement: " . $write_dbh->errstr;

my $i;
seek( $fh, 0, 0 );
$clusters_processed = 0;
my @clusters;
while (<$fh>) {
    $clusters_processed++;

    chomp;

    my @values = split( ',', $_ );

    # warn "Line: " . join ":", @values;

    my $repeatcount = 0;
    my $refcount    = 0;
    my $readcount   = 0;
    my $minpat      = 1000000;
    my $maxpat      = 0;
    my $range       = 0;

    # process each line
    # for statistics
    foreach my $val (@values) {

        $val =~ s/[\'\"]//g;

        # insert clusterlnk entry
        if ( $val <= 0 ) {

            $refcount++;

            $totalRefReps++;

        }
        else {

            $readcount++;

            $totalReadReps++;

        }

        $repeatcount++;
    }

    # execute ref and read pulls
    $sth->execute($clusters_processed)
        or die "Couldn't execute statement: " . $sth->errstr;
    $sth1->execute($clusters_processed)
        or die "Couldn't execute statement: " . $sth1->errstr;

    # store refs for later use
    open( my $RFILE, ">$TEMPDIR/refs_$DBSUFFIX.txt" ) or die $!;
    while ( my @data = $sth->fetchrow_array() ) {
        print $RFILE "-"
            . $data[0]
            . $data[6] . ","
            . $data[1] . ","
            . $data[2] . ","
            . $data[3] . ","
            . $data[4] . "\n";

        $minpat = min( $minpat, length( $data[4] ) );
        $maxpat = max( $maxpat, length( $data[4] ) );

        $i++;
    }
    close($RFILE);

    # print reads in blocks of BREAK_SIZE
    $i = 0;
    my $hcount = 0;
    while ( my @data = $sth1->fetchrow_array() ) {
        if ( ( $i % $BREAK_SIZE ) == 0 ) {
            $hcount++;

            print "@($clusters_processed\_$hcount):";
            print
                "\n**********************************************************************\n";

            # print refs each time
            open( my $RFILE, "<$TEMPDIR/refs_$DBSUFFIX.txt" ) or die $!;
            while (<$RFILE>) { print $_; }
            close($RFILE);
        }
        my $dna = nowhitespace( $data[1] );
        print $data[0]
            . $data[6] . ","
            . $data[2] . ","
            . $data[3] . ","
            . $dna . ","
            . $data[4] . "\n";

        $minpat = min( $minpat, length( $data[4] ) );
        $maxpat = max( $maxpat, length( $data[4] ) );

        $i++;
    }

    # do for 1st 10 clusters for now
    #if ($clusters_processed >= 20) { last; }

    # insert database records (cluster table)
    if ( $ENV{DEBUG} ) {
        my $numrefs  = $sth->rows;
        my $numreads = $sth1->rows;
        warn
            "Cluster $clusters_processed, numrefs: $numrefs, numreads: $numreads\n";
    }
    push @clusters,
        [ $clusters_processed, $minpat, $maxpat, $repeatcount, $refcount ];

    # $sth2->execute( $clusters_processed, $minpat, $maxpat, $repeatcount,
    #     $refcount )    # Execute the query
    #     or die "Couldn't execute statement: " . $sth2->errstr;
    if ( ( @clusters % $RECORDS_PER_INFILE_INSERT == 0 ) ) {
        my $cb   = gen_exec_array_cb( \@clusters );
        my $rows = vs_db_insert( $write_dbh, $sth2, $cb,
            "Error when inserting entries into our clusters table.\n" );
        if ($rows) {
            @clusters = ();
        }
        else {
            die
                "Something went wrong inserting, but somehow wasn't caught!\n";
        }
    }

    # $dbh->commit;

    # stats
    $range       = int( ( $maxpat / $minpat - 1.0 ) * 100 + .5 );
    $mostReps    = max( $mostReps, $repeatcount );
    $mostRefReps = max( $mostRefReps, $refcount );
    $maxRange    = max( $maxRange, $range );

}    # end of while loop

if (@clusters) {
    my $cb   = gen_exec_array_cb( \@clusters );
    my $rows = vs_db_insert( $write_dbh, $sth2, $cb,
        "Error when inserting entries into our clusters table.\n" );
    if ($rows) {
        @clusters = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

# enable old settings
$write_dbh->do("PRAGMA foreign_keys = ON");
$write_dbh->do("PRAGMA synchronous = ON");

$write_dbh->disconnect();
$read_dbh->disconnect();

# update the stats table
my %stats = (
    CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER  => $mostReps,
    CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER => $mostRefReps,
    CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER       => $maxRange,
    CLUST_NUMBER_OF_PROCLU_CLUSTERS                => $clusters_processed,
    CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS           => $totalRefReps,
    CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS          => $totalReadReps,
);
set_statistics( \%stats );

print STDERR
    "Processing complete -- processed $clusters_processed cluster(s)."
    . strftime( "\n\nend: %F %T\n\n", localtime );

