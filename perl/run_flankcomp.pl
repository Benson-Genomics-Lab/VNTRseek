#!/usr/bin/env perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use strict;
use warnings;
use Cwd;
use DBI;
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil
    qw(get_config get_dbh set_statistics gen_exec_array_cb vs_db_insert);


my $argc = @ARGV;
die "Usage: run_flankcomp.pl expects 4 arguments.\n"
    unless $argc >= 4;

my $curdir    = getcwd();
my $inputfile = $ARGV[0];
my $cnf       = $ARGV[1];
my $TEMPDIR   = $ARGV[2];
my $outfile   = $ARGV[3];

# get run config
my %run_conf = get_config("CONFIG", $cnf);
my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } )
    or die "Could not connect to database for reading: $DBI::errstr";
my $write_dbh = get_dbh()
    or die "Could not connect to database for writing: $DBI::errstr";

my $clusters_processed = 0;
my $totalRefReps       = 0;
my $totalReadReps      = 0;
my $sth;
my $sth1;
my $sth2;

my $mostReps    = 0;
my $mostRefReps = 0;
my $maxRange    = 0;

my $BREAK_SIZE = 4000;

$write_dbh->do("PRAGMA foreign_keys = OFF");
$write_dbh->do("PRAGMA synchronous = OFF");

# clear database cluster tables
$write_dbh->do( "DELETE FROM clusters" )
    or die "Couldn't do statement: " . $write_dbh->errstr;

$write_dbh->do( "DELETE FROM clusterlnk" )
    or die "Couldn't do statement: " . $write_dbh->errstr;


#############################################################################################

print "Inserting into clusterlnk table.\n";

open my $fh, "<$inputfile" or die $!;

$sth = $write_dbh->prepare(q{INSERT INTO clusterlnk VALUES (?, ?, ?, 0, 0)})
    or die "Couldn't prepare statement: " . $write_dbh->errstr;

# insert into clusterlnk
my $totalreps = 0;

# DBI->trace("3|SQL", "dbitrace.log"); # what is this wizardry?
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
            @cluster_links = ();
        }

    }

}    # end of while loop

# Finish insert
if (@cluster_links) {
    my $cb   = gen_exec_array_cb( \@cluster_links );
    my $rows = vs_db_insert( $write_dbh, $sth, $cb,
        "Error when inserting entries into our clusterlnk table.\n" );
    @cluster_links = ();
}
$sth->finish();

#############################################################################################

# Sync the DB so the next SELECTs work
$write_dbh->do("PRAGMA synchronous = FULL"); # KA: i'm confident this doesnt matter
$write_dbh->do("PRAGMA synchronous = OFF");

print "Printing DNA and inserting into cluster table.\n";

# now print dna and quals (also insert into cluster table)
$sth = $read_dbh->prepare(
    q{SELECT rid, flankleft, sequence, flankright, pattern, copynum, direction
    FROM refdb.fasta_ref_reps
    INNER JOIN clusterlnk ON rid=-repeatid
    WHERE clusterid = ?})
    or die "Couldn't prepare statement: " . $read_dbh->errstr;
$sth1 = $read_dbh->prepare(
    q{SELECT rid, dna, first, last, pattern, copynum, direction
    FROM fasta_reads
    INNER JOIN replnk ON fasta_reads.sid=replnk.sid
    INNER JOIN clusterlnk ON rid=repeatid
    WHERE clusterid = ?})
    or die "Couldn't prepare statement: " . $read_dbh->errstr;
$sth2 = $write_dbh->prepare(
    q{INSERT INTO clusters(cid, minpat, maxpat, repeatcount, refcount)
    VALUES(?,?,?,?,?)})
    or die "Couldn't prepare statement: " . $write_dbh->errstr;

my $i;
seek( $fh, 0, 0 );
$clusters_processed = 0;
my @clusters;
while (<$fh>) {
    $clusters_processed++;

    chomp;
    my @values = split( ',', $_ );

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
    open( my $RFILE, ">$TEMPDIR/refs.txt" ) or die $!;
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

    open( my $outfh, ">>$outfile" ) or die $!;
    while ( my @data = $sth1->fetchrow_array() ) {
        if ( ( $i % $BREAK_SIZE ) == 0 ) {
            $hcount++;

            print $outfh "@($clusters_processed\_$hcount):";
            print $outfh
                "\n**********************************************************************\n";

            # print refs each time
            open( my $RFILE, "<$TEMPDIR/refs.txt" ) or die $!;
            while (<$RFILE>) { print $outfh $_; }
            close($RFILE);
        }
        my $dna = nowhitespace( $data[1] );
        print $outfh $data[0]
            . $data[6] . ","
            . $data[2] . ","
            . $data[3] . ","
            . $dna . ","
            . $data[4] . "\n";

        $minpat = min( $minpat, length( $data[4] ) );
        $maxpat = max( $maxpat, length( $data[4] ) );

        $i++;
    }
    close($outfh);

    # insert database records (cluster table)
    if ( $ENV{DEBUG} ) {
        my $numrefs  = $sth->rows;
        my $numreads = $sth1->rows;
        warn "Cluster $clusters_processed, numrefs: $numrefs, numreads: $numreads\n";
    }
    push @clusters,
        [ $clusters_processed, $minpat, $maxpat, $repeatcount, $refcount ];

    # $sth2->execute( $clusters_processed, $minpat, $maxpat, $repeatcount,
    #     $refcount )    # Execute the query
    #     or die "Couldn't execute statement: " . $sth2->errstr;
    if ( ( @clusters % $RECORDS_PER_INFILE_INSERT == 0 ) ) {
        my $cb   = gen_exec_array_cb( \@clusters );
        my $rows = vs_db_insert( $write_dbh, $sth2, $cb,
            "Error when inserting entries into clusters table.\n" );
        @clusters = ();
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
        "Error when inserting entries into clusters table.\n" );
    @clusters = ();
}
$sth->finish();
$sth1->finish();
$sth2->finish();

# enable old settings
$write_dbh->do("PRAGMA foreign_keys = ON");
$write_dbh->do("PRAGMA synchronous = ON");
$write_dbh->disconnect();

$read_dbh->disconnect();

# update the stats table
set_statistics({
    CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER  => $mostReps,
    CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER => $mostRefReps,
    CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER       => $maxRange,
    CLUST_NUMBER_OF_PROCLU_CLUSTERS                => $clusters_processed,
    CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS           => $totalRefReps,
    CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS          => $totalReadReps,
});

print "Processing complete -- processed $clusters_processed cluster(s).";

1;

sub nowhitespace {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}
