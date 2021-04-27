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
use Try::Tiny;

use lib "$FindBin::RealBin/lib";

use vutil
    qw(get_config get_dbh set_statistics get_trunc_query gen_exec_array_cb vs_db_insert);

my $updatedClustersCount = 0;
my $updatedRefsCount     = 0;

sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $curdir = getcwd;

my $argc = @ARGV;
if ( $argc < 5 ) {
    die
        "Usage: run_rankflankmap.pl inputfile  mapdir tmpdir dbsuffix run_dir\n";
}

my $inputfile = $ARGV[0];
my $mapdir    = $ARGV[1];
my $tmp       = $ARGV[2];
my $DBSUFFIX  = $ARGV[3];
my $run_dir   = $ARGV[4];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $DBSUFFIX, $run_dir );

my $clusters_processed = 0;

my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth1;
my $map_insert_sth;
my $rankflank_insert_sth;

#goto AAA;
# disable indices
$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");
$dbh->do("PRAGMA journal_mode = TRUNCATE");

# clear map
$dbh->begin_work;
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "map" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

# clear rankflank
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "rankflank" ) )
    or die "Couldn't do statement: " . $dbh->errstr;
$dbh->commit;

$map_insert_sth = $dbh->prepare(
    qq{INSERT INTO map (refid, readid, reserved, reserved2)
    VALUES (?, ?, 0, 0)}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$rankflank_insert_sth
    = $dbh->prepare(qq{INSERT INTO rankflank VALUES (?, ?, ?, ?)})
    or die "Couldn't prepare statement: " . $dbh->errstr;

open FILE, "<$inputfile" or die "error opening for reading '$inputfile': $!";

opendir( DIR, $mapdir );
my @allfiles = readdir(DIR);
closedir(DIR);

my $j = 0;
my $k = 0;
my $upload;
my $uploadedrank = 0;
my $uploadedmap  = 0;
my ( @map_rows, @rankflank_rows );

#foreach my $file (@files) {
foreach my $file (@allfiles) {

    if ( index( $file, "_map" ) ) {

        $clusters_processed++;

        my $mfile = "$mapdir/$file";
        open MFILE, "<$mfile" or die $!;

        my %READVECTOR = ();
        while (<MFILE>) {
            chomp;
            my @mfields = split( '\t', $_ );
            my $msize   = scalar @mfields;
            if ( $msize >= 8 ) {
                my $readid = $mfields[0];

                #print "\n\n".$_."\n\n";
                #exit(1);

                my @rfields = split( ',', $mfields[7] );

                #my @temparray = ();

                my $bestscore = 0;
                my $bestref   = "";

                foreach my $refstr (@rfields) {
                    if ( $refstr =~ /^-?(\d+):(\d+):(\d+)/ ) {

                        #print "\n".$1." ".$2." ".$3."\n";
                        #exit(1);

                        # calculate score
                        my $score;

                        if ( ( $mfields[2] + $mfields[3] )
                            == 0
                            ) # if no flanks, it will only be marked best if nothing else is available
                        {
                            $score = 0;
                        }
                        else

                   #{ # asked to add by Dr. Benson to balance out small flanks

                         #  my $A = ($2 + $3) / ( $mfields[2] + $mfields[3] );
                         #  my $B = 1 / ( $mfields[2] + $mfields[3]);
                         #  $score = 1 - max($A,$B);
                         #}

                        {
                            $score = 1 - (
                                ( $2 + $3 ) / ( $mfields[2] + $mfields[3] ) );
                        }

               # filter to remove all flank scores below .9, added Nov 5, 2012
                        if ( $score >= 0.90 ) {

                            if ( $score > $bestscore ) {
                                $bestref   = $1;
                                $bestscore = $score;
                            }
                            elsif ( $score == $bestscore ) {
                                if ( $bestref eq "" ) { $bestref = $1; }
                                else { $bestref .= ( "," . $1 ); }
                            }
                        }

       # create a map entry in database (added 11/19/2010)
       #$map_insert_sth->execute($1,$readid)             # Execute the query
       #      or die "Couldn't execute statement: " . $map_insert_sth->errstr;
                        $k++;
                        push @map_rows, [ $1, $readid ];
                        if ( @map_rows % $RECORDS_PER_INFILE_INSERT == 0 ) {
                            my $cb = gen_exec_array_cb( \@map_rows );
                            my $rows
                                = vs_db_insert( $dbh, $map_insert_sth, $cb,
                                "Error inserting into map table." );
                            if ($rows) {
                                $uploadedmap += $rows;
                                @map_rows = ();
                            }
                            else {
                                die
                                    "Something went wrong inserting, but somehow wasn't caught!\n";
                            }
                        }
                    }
                }

                # insert the rankflank
                if ( $bestref ne "" ) {
                    my @ranks = split( ',', $bestref );
                    my $ties  = scalar(@ranks) - 1;
                    foreach my $rstr (@ranks) {

#$rankflank_insert_sth->execute($readid,$rstr,$bestscore)             # Execute the query
#      or die "Couldn't execute statement: " . $rankflank_insert_sth->errstr;
                        $j++;

                        push @rankflank_rows,
                            [ $rstr, $readid, $bestscore, $ties ];
                        if ( @rankflank_rows % $RECORDS_PER_INFILE_INSERT
                            == 0 )
                        {
                            my $cb = gen_exec_array_cb( \@rankflank_rows );
                            my $rows
                                = vs_db_insert( $dbh, $rankflank_insert_sth,
                                $cb,
                                "Error inserting into rankflank table." );
                            if ($rows) {
                                $uploadedrank += $rows;
                                @rankflank_rows = ();
                            }
                            else {
                                die
                                    "Something went wrong inserting, but somehow wasn't caught!\n";
                            }
                        }

                    }
                }

                #$READVECTOR{$readid} = @temparray;

                #print "\n";
            }
        }

        # load the files and remove the temp files

        ( $ENV{DEBUG} ) && warn "processed: $clusters_processed\n";

    }    # end of if (index($file,"_map")) {

}    # end of foreach @files
close(FILE);

if (@map_rows) {
    my $cb   = gen_exec_array_cb( \@map_rows );
    my $rows = vs_db_insert( $dbh, $map_insert_sth, $cb,
        "Error inserting into map table." );
    if ($rows) {
        $uploadedmap += $rows;
        @map_rows = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

if ( $uploadedmap != $k ) {
    die
        "Uploaded number of map entries($uploadedmap) not equal to the number of uploaded counter ($k), aborting!\n";
}

if (@rankflank_rows) {
    my $cb   = gen_exec_array_cb( \@rankflank_rows );
    my $rows = vs_db_insert( $dbh, $rankflank_insert_sth, $cb,
        "Error inserting into rankflank table." );
    if ($rows) {
        $uploadedrank += $rows;
        @rankflank_rows = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

if ( $uploadedrank != $j ) {
    die
        "Uploaded number of rankflank entries($uploadedrank) not equal to the number of uploaded counter ($j), aborting!\n";
}
AAA:

# create temp table for deletions
$dbh->do(
    'CREATE TEMPORARY TABLE ranktemp (
    `refid` integer NOT NULL,
    `readid` integer NOT NULL,
    PRIMARY KEY (`refid`, `readid`))'
) or die "Couldn't do statement: " . $dbh->errstr;

$sth1 = $dbh->prepare(qq{INSERT INTO ranktemp VALUES (?, ?)});
print STDERR
    "Prunning (keep best ref for each read) from rankflank table...\n";
my $query
    = "Select refid,readid,sid,score FROM rankflank INNER JOIN replnk ON rankflank.readid=replnk.rid ORDER BY readid, score;";
$sth = $dbh->prepare($query);
$sth->execute();
my $i        = 0;
my $count    = 0;
my $oldseq   = -1;
my $oldref   = -1;
my $oldread  = -1;
my $oldscore = -1.0;

my @rows;
while ( my @data = $sth->fetchrow_array() ) {
    if ( $data[1] == $oldread && $data[3] != $oldscore ) {

        # delete old one
        push @rows, [ $oldref, $oldread ];
        if ( @rows % $RECORDS_PER_INFILE_INSERT == 0 ) {
            my $cb  = gen_exec_array_cb( \@rows );
            my $num = vs_db_insert( $dbh, $sth1, $cb,
                "Error inserting into temporary rank table." );
            if ($num) {
                $count += $num;
                @rows = ();
            }
            else {
                die
                    "Something went wrong inserting, but somehow wasn't caught!\n";
            }
        }
    }
    $oldref   = $data[0];
    $oldread  = $data[1];
    $oldseq   = $data[2];
    $oldscore = $data[3];
    $i++;
}

if (@rows) {
    my $cb  = gen_exec_array_cb( \@rows );
    my $num = vs_db_insert( $dbh, $sth1, $cb,
        "Error inserting into temporary rank table." );
    if ($num) {
        $count += $num;
        @rows = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

print STDERR "Prunning complete. Pruned $count rankflank records.\n";

print STDERR "Prunning all (one TR/same read) rankflank table...\n";
$query
    = "Select refid,readid,sid,score FROM rankflank INNER JOIN replnk ON rankflank.readid=replnk.rid ORDER BY refid, sid, score, readid;"
    ; # readid added for tie resolution to keep rank and rankflank entries more in sync
$sth = $dbh->prepare($query);
$sth->execute();
$i       = 0;
$count   = 0;
$oldseq  = -1;
$oldref  = -1;
$oldread = -1;

while ( my @data = $sth->fetchrow_array() ) {

    if ( $data[0] == $oldref && $data[2] == $oldseq ) {

        push @rows, [ $oldref, $oldread ];
        if ( @rows % $RECORDS_PER_INFILE_INSERT == 0 ) {
            my $cb  = gen_exec_array_cb( \@rows );
            my $num = vs_db_insert( $dbh, $sth1, $cb,
                "Error inserting into temporary rank table." );
            if ($num) {
                $count += $num;
                @rows = ();
            }
            else {
                die
                    "Something went wrong inserting, but somehow wasn't caught!\n";
            }
        }
    }
    $oldref  = $data[0];
    $oldread = $data[1];
    $oldseq  = $data[2];
    $i++;
}

if (@rows) {
    my $cb  = gen_exec_array_cb( \@rows );
    my $num = vs_db_insert( $dbh, $sth1, $cb,
        "Error inserting into temporary rank table." );
    if ($num) {
        $count += $num;
        @rows = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}

print STDERR "Prunning complete. Pruned $count rankflank records.\n";

# delete from rankflank based on temptable entries
my $delfromtable = 0;
$query = qq{
    DELETE FROM rankflank
    WHERE EXISTS (
        SELECT * FROM ranktemp t2
        WHERE rankflank.refid = t2.refid
            AND rankflank.readid = t2.readid
    )
};
$dbh->begin_work;
$delfromtable = $dbh->do($query);
$dbh->commit;

# $sth->finish;

# set old settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");

if ( $delfromtable != $count ) {
    die
        "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($count), aborting!";
}

$dbh->disconnect();
set_statistics(
    {   RANKFLANK_EDGES_INSERTED  => $j,
        RANKFLANK_REMOVED_SAMEREF => $count,
        RANKFLANK_REMOVED_SAMESEQ => $count
    }
);

print STDERR "\n\n";

print STDERR
    "Processing complete -- processed $clusters_processed cluster(s). Deleted from rankflank using temptable: $delfromtable\n";

warn strftime( "\n\nend: %F %T\n\n", localtime );

1;

