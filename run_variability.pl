#!/usr/bin/env perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use strict;
use warnings;
use List::Util qw[min max];
use Cwd;
use DBI;
use POSIX qw(strftime);
use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil
    qw(get_config get_dbh set_statistics gen_exec_array_cb vs_db_insert get_trunc_query);

if ( $ENV{DEBUG} ) {
    use Data::Dumper;
}

sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

my $updatedClustersCount = 0;
my $updatedRefsCount     = 0;

my $curdir = getcwd;

my $argc = @ARGV;
if ( $argc < 6 ) {
    die
        "Usage: run_variability.pl inputfile  mapdir dbname run_dir minflank tempdir\n";
}

my $inputfile          = $ARGV[0];
my $mapdir             = $ARGV[1];
my $DBSUFFIX           = $ARGV[2];
my $run_dir            = $ARGV[3];
my $MIN_FLANK_REQUIRED = $ARGV[4];
my $TEMPDIR            = $ARGV[5];

warn strftime( "\n\nstart: %F %T\n\n", localtime );
my %run_conf = get_config( $DBSUFFIX, $run_dir );
my $dbh      = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

my ( $sth, $sth1, $sth6, $sth7, $sth8, $query );
my $TEMPFILE;
my $TEMP_CLNK;

# change settings to speedup updates and inserts
$dbh->do("PRAGMA foreign_keys = OFF");

$dbh->begin_work;

# update reserved field on entire table
$dbh->do('UPDATE clusterlnk SET reserved=0,reserved2=0;')
    or die "Couldn't do statement: " . $dbh->errstr;
$dbh->do('UPDATE map SET reserved=0,reserved2=0;')
    or die "Couldn't do statement: " . $dbh->errstr;
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "vntr_support" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->do(
    q{CREATE TEMPORARY TABLE mapr (
  `refid` INT(11) NOT NULL,
  `readid` INT(11) NOT NULL,
  PRIMARY KEY (refid,readid))}
);

# create temp table for clusterlnk table updates and update clusterlnk table
$dbh->do(
    q{CREATE TEMPORARY TABLE ctrlnk (
    `clusterid` INT(11) NOT NULL,
    `repeatid` INT(11) NOT NULL,
    `change` INT(11) NOT NULL,
    PRIMARY KEY (clusterid,repeatid)
    )}
);
$dbh->commit;

# prepare statments
my $read_dbh = get_dbh( { userefdb => 1, readonly => 1 } );
$sth = $read_dbh->prepare(
    q{SELECT rid,flankleft,sequence,flankright,pattern,copynum,(lastindex-firstindex+1) AS arlen
  FROM refdb.fasta_ref_reps
  WHERE rid = ?}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth1 = $read_dbh->prepare(
    q{SELECT rid, dna, first, last, pattern, copynum
  FROM fasta_reads INNER JOIN replnk ON fasta_reads.sid=replnk.sid
  WHERE rid = ?}
) or die "Couldn't prepare statement: " . $dbh->errstr;

$sth6 = $dbh->prepare(q{INSERT INTO ctrlnk VALUES(?, ?, ?)})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth7 = $dbh->prepare(q{INSERT INTO mapr VALUES(?, ?)})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth8 = $read_dbh->prepare(
    q{SELECT map.refid, map.readid from clusterlnk
  INNER JOIN map ON map.refid=-clusterlnk.repeatid
  WHERE map.bbb=1 AND clusterlnk.clusterid=?
  ORDER BY map.readid,map.refid}
) or die "Couldn't prepare statement: " . $dbh->errstr;

my ( %VNTR_REF, %VNTR_COPIES, %VNTR_COPIESFLOAT, %VNTR_SUPPORT,
    %VNTR_SAMEASREF, %VNTR_REPRESENTATIVE );

my $SUPPORT_INCREMENTED = 0;

my $processed    = 0;
my $cl_processed = 0;

open my $fh, "<$inputfile" or die $!;

my $clusters_processed = 0;
while (<$fh>) {
    $clusters_processed++;

    my @values = split( ',', $_ );

    my %REFCOPIES   = ();
    my %ASKLENGTH   = ();
    my %READCOPIES  = ();
    my $repeatcount = 0;
    my $refcount    = 0;
    my $readcount   = 0;

    my $readlen;
    my $first;
    my $last;

    foreach my $val (@values) {

        my $dir = q{'};
        if ( $val =~ m/([\'\"])/ ) { $dir = $1; }

        chomp($val);

        $val =~ s/[\'\"]//g;

        $repeatcount++;

        # go though all refs for each read
        if ( $val <= 0 ) {

            $refcount++;

            #$val = ($val<0)? -$val : $val; # ids are positive in database

            $sth->execute( -$val )    # Execute the query
                or die "Couldn't execute statement: " . $sth->errstr;

            my @data = $sth->fetchrow_array();

            if ( $sth->rows == 0 ) {
                print STDERR
                    "No record in database for entry `$val'. Aborting!\n\n";
                exit(1);
            }

            $ASKLENGTH{$val}
                = int( $data[6] ); # remember length of the ref array sequence
            $REFCOPIES{$val} = $data[5];

        }
        else {

            $readcount++;

            $sth1->execute($val)    # Execute the query
                or die "Couldn't execute statement: " . $sth1->errstr;

            my @data = $sth1->fetchrow_array();

            if ( $sth1->rows == 0 ) {
                print STDERR
                    "No record in database for entry `$val'. Aborting!\n";
                exit(1);
            }

            $readlen = length( nowhitespace( $data[1] ) );

            $ASKLENGTH{$val} = $readlen;    # remember length of the read

            $first = $data[2];
            $last  = $data[3];
            if (   ( $first - 1 ) >= $MIN_FLANK_REQUIRED
                && ( $readlen - $last ) >= $MIN_FLANK_REQUIRED )
            {

                $READCOPIES{$val} = $data[5];
            }

        }

    }

    # get map data
    my %READVECTOR = ();
    my %REFHASH    = ();
    my $readidold  = -1;
    $sth8->execute($clusters_processed)    # Execute the query
        or die "Couldn't execute statement: " . $sth8->errstr;
    while ( my @data = $sth8->fetchrow_array() ) {
        my $readid;
        my $refid;
        my $readlen;

        $refid  = -1 * int( $data[0] );
        $readid = int( $data[1] );

        $REFHASH{$refid} = 1;

        warn "\n\n$refid | $readid \n";

        if ( $readid != $readidold ) {
            $READVECTOR{$readid} = ();
        }

        $readlen = $ASKLENGTH{$readid};

        if ( !exists $ASKLENGTH{$readid} ) {
            warn "Entry for $readid does not exist!\n";
            exit(1);
        }
        if ( !exists $ASKLENGTH{$refid} ) {
            warn "Entry for $refid does not exist!\n";
            exit(1);
        }

        #if ($ASKLENGTH{$refid} <= $readlen) {
        push( @{ $READVECTOR{$readid} }, $refid );

        #}

        $readidold = $readid;
    }
    $sth8->finish;
    print STDERR "\nprocessed: $clusters_processed\n\n";

    # do for 1st 10 clusters for now
    # if ($clusters_processed >= 20) { last; }

    # insert database records (cluster table)
    #print "VYN:  ";
    #print  VNTR_YES_NO(\%REFCOPIES,\%READCOPIES,\%READVECTOR);

    my $vYes
        = VNTR_YES_NO( \$TEMPFILE, \$TEMP_CLNK, \%REFCOPIES, \%READCOPIES,
        \%READVECTOR, \%REFHASH, $clusters_processed );

    $updatedClustersCount += $vYes;

}

warn "VNTR_SUPPORT hash:\n", Dumper( \%VNTR_SUPPORT ), "\n"
    if ( $ENV{DEBUG} );

close($fh);
$sth->finish;
$sth1->finish;

# $dbh->begin_work;
my $updCLNKfromfile = $dbh->do(
    q{UPDATE clusterlnk SET reserved=(
    SELECT change FROM ctrlnk t2
    WHERE clusterlnk.clusterid = t2.clusterid
        AND clusterlnk.repeatid=t2.repeatid
    )
    WHERE EXISTS (SELECT * FROM ctrlnk t2
    WHERE clusterlnk.clusterid = t2.clusterid
        AND clusterlnk.repeatid=t2.repeatid)}
);

# update map based on temp table
my $updfromtable = $dbh->do(
    q{UPDATE map SET reserved=1
    WHERE EXISTS (
    SELECT * FROM mapr t2
    WHERE map.refid = t2.refid
        AND map.readid = t2.readid)}
);

# $dbh->commit;

# write SUPPORT info to temp files to be loaded into vntr_support

# ($updfromtable) = $dbh->selectrow_array(q{SELECT COUNT(*) FROM map WHERE reserved=1});
$query = q{INSERT INTO vntr_support VALUES(?, ?, ?, ?, ?, ?)};
$sth   = $dbh->prepare($query);

my $supcounter = 0;
my $supInsert  = 0;
my @vntrs_supported;
foreach my $key ( keys %VNTR_REF ) {
    $supcounter++;
    push @vntrs_supported,
        [
        $VNTR_REF{$key},
        $VNTR_COPIES{$key},
        $VNTR_SAMEASREF{$key},
        $VNTR_SUPPORT{$key},
        $VNTR_COPIESFLOAT{$key},
        ( exists $VNTR_REPRESENTATIVE{$key} )
        ? $VNTR_REPRESENTATIVE{$key}
        : undef
        ];
    if ( @vntrs_supported % $RECORDS_PER_INFILE_INSERT == 0 ) {
        my $cb = gen_exec_array_cb( \@vntrs_supported );
        my $rows
            = vs_db_insert( $dbh, $sth, $cb,
            "Error inserting into vntr_support table." );
        if ($rows) {
            $supInsert += $rows;
            @vntrs_supported = ();
        }
        else {
            die
                "Something went wrong inserting, but somehow wasn't caught!\n";
        }
    }
}

if (@vntrs_supported) {
    my $cb = gen_exec_array_cb( \@vntrs_supported );
    my $rows
        = vs_db_insert( $dbh, $sth, $cb,
        "Error inserting into vntr_support table." );
    if ($rows) {
        $supInsert += $rows;
        @vntrs_supported = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

$dbh->begin_work;
$dbh->do(
    q{CREATE TEMPORARY TABLE ctr (
    `clusterid` INT(11) NOT NULL PRIMARY KEY,
    `varbl` INT(11) NOT NULL DEFAULT 0
    )}
);

my $InsClusToFile = $dbh->do(
    q{INSERT INTO ctr SELECT clusterid, count(*) as vrefs
    FROM clusterlnk
    WHERE reserved>0
    GROUP by clusterid
    ORDER BY vrefs DESC}
);

my $UpdClusFromFile = $dbh->do(
    q{UPDATE OR IGNORE clusters
    SET variability=(SELECT varbl FROM ctr t2
    WHERE clusters.cid == t2.clusterid)}
);
$dbh->commit;

# set old settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->disconnect();

if ( $supInsert != $supcounter ) {
    die
        "Inserted number of vntr_support entries($supInsert) not equal to the number of inserted counter ($supcounter), aborting!";
}
if ( $updfromtable != $processed ) {
    die
        "Updated number of map entries($updfromtable) not equal to the number of inserted counter ($processed), aborting!";
}
if ( $updCLNKfromfile != $cl_processed ) {
    die
        "Updated number of cluster entries($updCLNKfromfile) not equal to the number of inserted counter ($cl_processed), aborting!";
}
if ( $UpdClusFromFile != $InsClusToFile ) {
    die
        "Updated number of clusterlnk entries($UpdClusFromFile) not equal to the number of inserted counter ($InsClusToFile), aborting!";
}

print STDERR "\n\n";

print STDERR
    "Processing complete -- processed $clusters_processed cluster(s), support entries created = $supInsert.\n";

set_statistics(
    {   CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR     => $updCLNKfromfile,
        CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR => $updatedClustersCount
    }
);

warn strftime( "\n\nend: %F %T\n\n", localtime );

1;

# this function adds support to reference copynumber count
sub add_support {

    my $refid          = $_[0];
    my $sameasref      = $_[1];
    my $copies         = $_[2];
    my $copiesfloat    = $_[3];
    my $representative = $_[4];

    my $key = $refid . "_" . $copies;

    if ( exists $VNTR_SUPPORT{$key} ) {

        $VNTR_SUPPORT{$key}++;

    }
    else {

        $VNTR_REF{$key}            = $refid;
        $VNTR_SAMEASREF{$key}      = $sameasref;
        $VNTR_COPIES{$key}         = $copies;
        $VNTR_COPIESFLOAT{$key}    = $copiesfloat;
        $VNTR_REPRESENTATIVE{$key} = $representative;
        $VNTR_SUPPORT{$key}        = 1;
    }

    $SUPPORT_INCREMENTED++;

    return 0;
}

# shows that there is 0 support, ignores if already an entry for this reference
sub add_zero_support {

    my $refid       = $_[0];
    my $sameasref   = $_[1];
    my $copies      = $_[2];
    my $copiesfloat = $_[3];

    my $key = $refid . "_" . $copies;

    if ( !exists $VNTR_SUPPORT{$key} ) {

        $VNTR_REF{$key}         = $refid;
        $VNTR_SAMEASREF{$key}   = $sameasref;
        $VNTR_COPIES{$key}      = $copies;
        $VNTR_COPIESFLOAT{$key} = $copiesfloat;
        $VNTR_SUPPORT{$key}     = 0;
    }

    return 0;
}

# this is a function that detects if reads have (somewhat) different number of copy numbers then reference(s)
# TODO Rewrite signature for this function: we'll no longer use mapr_fh and clnk_fh
sub VNTR_YES_NO {

    my $mapr_fh  = ${ shift() };
    my $clnk_fh  = ${ shift() };
    my %refs     = %{ shift() };
    my %reads    = %{ shift() };
    my %readhash = %{ shift() };
    my %newrefs  = %{ shift() }
        ; # this is to make sure we have an etnry in vntr_support for all BBB refs
    my $clusterid = shift;

    my $varyes = 0;
    my $change = 0;

    my (%REF_UPDATED) = ();

    #print "\nVNTR_YES_NO:";

    # for each read in map file
    # warn Dumper( \%readhash );
    my @rows;
    while ( my ( $key, @temp ) = each(%readhash) ) {

        my $valread;
        my $valref;

        if ( exists $reads{$key} ) {

            $valread = $reads{$key};

            #print STDERR "\n$key($valread) =>";

    # for each references listed in map file as being associated with the read
    #foreach my $val (@{$temp[0]}) {
            foreach my $val ( @{ $readhash{$key} } ) {

                if ( exists $refs{$val} ) {
                    $valref = $refs{$val};
                }
                else {
                    print "Undefined reference `$val`. Aborting!\n";
                    exit(1);
                }

                #print STDERR " ".$val."($valref)";

              # print "$valRead - $valRef\n";
              #if (($valread >= ($valref-.5)) && ($valread <= ($valref+.5))) {
                if (   ( $valread > ( $valref - .8 ) )
                    && ( $valread < ( $valref + .8 ) ) )
                {

                    # add support
                    #add_support($val, 1, int(abs($valref) + 0.5), $valref );
                    add_support( $val, 1, 0, $valref, $key );

                }
                else {

                    $varyes = 1;

                    # add support
                    #add_support($val, 0, int(abs($valread) + 0.5) , 0.0);

                    my $fchange;

                    if ( $valread > $valref ) {
                        $fchange = 1 + int( $valread - ( $valref + .8 ) );
                    }    # truncation
                    else {
                        $fchange = -1 - int( ( $valref - .8 ) - $valread );
                    }    # truncation

                    add_support( $val, 0, $fchange, 0.0, $key );

            # this ignored if there is already an entry for this ref
            # this would be incremented to 0 if there is a read found later on
            # add_zero_support($val, 1, int(abs($valref) + 0.5) );

                    #print "\t$clusterid: $refid\n";

                    $change = int( abs( $valread - $valref ) + 0.5 );

                    if ( !exists $REF_UPDATED{$val}
                        || $change > $REF_UPDATED{$val} )
                    {

#$sth3->execute($change,$clusterid,$val)             # set the reserved field on reference
#   or die "Couldn't execute statement: " . $sth3->errstr;

                        $REF_UPDATED{$val} = $change;
                    }

#$sth7->execute(-$val,$key)             # set the reserved field on map read, (added 11/19/2010)
#    or die "Couldn't execute statement: " . $sth7->errstr;

                    if ( $ENV{DEBUG} ) {
                        warn "mapr insert: " . -$val . ", $key\n";
                    }
                    push @rows, [ -$val, $key ];
                    $processed++;

                    if ( @rows % $RECORDS_PER_INFILE_INSERT == 0 ) {
                        my $cb   = gen_exec_array_cb( \@rows );
                        my $rows = vs_db_insert( $dbh, $sth7, $cb,
                            "Error when inserting entries into our mapr table.\n"
                        );
                        if ($rows) {
                            @rows = ();
                        }
                        else {
                            die
                                "Something went wrong inserting, but somehow wasn't caught!\n";
                        }
                    }
                }

            }

            #print "\n";

        }    # end of exists

    }    # end of loop

    # TODO Optimize for SQLite: this function is called for every cluster,
    # might be better to preserve array for accumulating rows across calls
    if (@rows) {
        my $cb   = gen_exec_array_cb( \@rows );
        my $rows = vs_db_insert( $dbh, $sth7, $cb,
            "Error when inserting entries into our mapr table.\n" );
        if ($rows) {
            @rows = ();
        }
        else {
            die
                "Something went wrong inserting, but somehow wasn't caught!\n";
        }
    }

    # create self support for refs
    foreach my $val ( keys %newrefs ) {
        my $valref = $refs{$val};

        #add_zero_support($val, 1, int(abs($valref) + 0.5), $valref );
        add_zero_support( $val, 1, 0, $valref );
    }

    # write to clustmp file (to update clusterlnk entries)
    @rows = ();
    foreach my $val ( keys %REF_UPDATED ) {
        $cl_processed++;
        $change = $REF_UPDATED{$val};
        push @rows, [ $clusterid, $val, $change ];
        if ( $ENV{DEBUG} ) {
            warn "ctrlnk insert: $clusterid, $val, $change\n";
        }
    }

    # TODO: See TODO above for optimizing for SQLite
    if (@rows) {
        my $cb   = gen_exec_array_cb( \@rows );
        my $rows = vs_db_insert( $dbh, $sth6, $cb,
            "Error when inserting entries into our ctrlnk table.\n" );
        if ($rows) {
            @rows = ();
        }
        else {
            die
                "Something went wrong inserting, but somehow wasn't caught!\n";
        }
    }

    # print "\nNOT VARIABLE\n\n";
    return $varyes;

}    # end of func

