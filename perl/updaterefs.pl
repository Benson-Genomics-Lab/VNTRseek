#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use POSIX "strftime";
use File::Basename;
use Try::Tiny;
use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil
    qw(get_config get_dbh set_statistics get_statistics gen_exec_array_cb vs_db_insert);

=head1

updaterefs.pl - Calculates VNTRs from database and writes
final reports and VCF files

=cut

# Arguments
my $argc = @ARGV;
die "Usage: updaterefs.pl expects 4 arguments.\n"
    unless $argc >= 4;

my $curdir        = getcwd();
my $RUN_NAME      = $ARGV[0]; # dbname
my $cnf           = $ARGV[1]; # config_file
my $result_prefix = $ARGV[2]; # base of output file names
my $VERSION       = $ARGV[3]; # VERSION

my %run_conf = get_config("CONFIG", $cnf);
my ($MIN_SUPPORT_REQUIRED)
    = @run_conf{qw(MIN_SUPPORT_REQUIRED)};

my $RECORDS_PER_INFILE_INSERT = 100000;
############################ Procedures ########################################

=head1 FUNCTIONS

=head2 RC( sequence )

Returns the reverse complement of a DNA sequence

=cut

sub RC {

    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ACGT/TGCA/;
    return $seq;
}

=head2 make_vcf_rec( $tr_hash )

Takes a hash representing a supported TR record, and
returns a VCF record as a string.

=cut

sub make_vcf_rec {
    my $supported_tr = shift;

    my $qual = ".";
    my $filter = ( $supported_tr->{is_singleton} == 1 ) ? "PASS" : "SC";
    my $info = sprintf(
        "RC=%.2lf;RPL=%d;RAL=%d;RCP=%s;",
        $supported_tr->{copynum}, length( $supported_tr->{consenuspat} ),
        $supported_tr->{arlen},   $supported_tr->{consenuspat}
    );
    my $format = "GT:SP:CGL";

    my $vcf_rec = join(
        "\t",
        $supported_tr->{head},
        ( $supported_tr->{firstindex} - 1 ),
        "td" . $supported_tr->{rid},
        $supported_tr->{sequence},
        $supported_tr->{alt},
        $qual, $filter, $info, $format,
        join( ":",
            $supported_tr->{gt_string}, $supported_tr->{support},
            $supported_tr->{cgl} )
    ) . "\n";

    return $vcf_rec;
}

####################################
# Takes a boolean as an argument. If the boolean is anything perl considers false,
# then this function will only produce a VCF file for supported VNTRs. If true, all
# supported TRs are included in the file.

sub print_vcf {

    # Get needed stats
    my @stats = qw(PARAM_TRF
        FILE_REFERENCE_SEQ
        FILE_REFERENCE_LEB
        NUMBER_REF_TRS
        FOLDER_FASTA
        FOLDER_PROFILES
        NUMBER_READS
        NUMBER_TRS_IN_READS);
    my $stat_hash = get_statistics(@stats);
    my $dbh       = get_dbh( { readonly => 1, userefdb => 1 } );

    # Get total number of TRs supported
    my $numsup_sth = $dbh->prepare(q{
        SELECT count(distinct vntr_support.refid)
        FROM vntr_support
        WHERE support >= $MIN_SUPPORT_REQUIRED})
        or die "Couldn't prepare statement: " . $dbh->errstr;

    my $numsup = 0;
    $numsup_sth->execute() or die "Cannot execute: " . $numsup_sth->errstr();
    $numsup_sth->bind_columns( \$numsup );
    $numsup_sth->fetch();
    if ( !defined $numsup ) {
        die "Error getting number of supported TRs: " . $dbh->errstr;
    }
    $numsup_sth->finish();

    # Get number of VNTRs supported
    my $numvntrs_sth = $dbh->prepare(q{
        SELECT count(*)
        FROM fasta_ref_reps
        WHERE support_vntr > 0})
        or die "Couldn't prepare statement: " . $dbh->errstr;

    my $numvntrs = 0;
    $numvntrs_sth->execute()
        or die "Cannot execute: " . $numvntrs_sth->errstr();
    $numvntrs_sth->bind_columns( \$numvntrs );
    $numvntrs_sth->fetch;
    $numvntrs_sth->finish;
    if ( !defined $numsup ) {
        die "Error getting number of supported VNTRs: " . $dbh->errstr;
    }

    # $update_spanN_sth->finish;
    my $spanN_fn = "${result_prefix}.span${MIN_SUPPORT_REQUIRED}.vcf";
    open my $spanN_vcffile, ">", $spanN_fn
        or die "Can't open for writing $spanN_fn\n\n";

    my $allwithsupport_fn
        = "${result_prefix}.allwithsupport.span${MIN_SUPPORT_REQUIRED}.vcf";
    open my $allwithsupport_vcffile, ">", $allwithsupport_fn
        or die "Can't open for writing $allwithsupport_fn\n\n";

    my $vcf_header
        = "##fileformat=VCFv4.2\n"
        . strftime( "##fileDate=\"%Y%m%d\"\n", localtime )
        . qq[##source=Vntrseek ver. $VERSION
##TRFParameters=$stat_hash->{PARAM_TRF}
##referenceseq=$stat_hash->{FILE_REFERENCE_SEQ}
##referenceprofile=$stat_hash->{FILE_REFERENCE_LEB}
##ploidy=$run_conf{PLOIDY}
##numrefTRs=$stat_hash->{NUMBER_REF_TRS}
##readseqfolder=$stat_hash->{FOLDER_FASTA}
##readprofilefolder=$stat_hash->{FOLDER_PROFILES}
##numreads=$stat_hash->{NUMBER_READS}
##numreadTRs=$stat_hash->{NUMBER_TRS_IN_READS}
##numVNTRs=$numvntrs
##numTRsWithSupport=$numsup
##database=${RUN_NAME}_rl$run_conf{READ_LENGTH}
##INFO=<ID=RC,Number=1,Type=Float,Description="Reference Copies">
##INFO=<ID=RPL,Number=1,Type=Integer,Description="Reference Pattern Length">
##INFO=<ID=RAL,Number=1,Type=Integer,Description="Reference Tandem Array Length">
##INFO=<ID=RCP,Number=1,Type=String,Description="Reference Consensus Pattern">
##INFO=<ID=ALGNURL,Number=1,Type=String,Description="Alignment URL">
##FILTER=<ID=SC,Description="Reference is Singleton">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SP,Number=R,Type=Integer,Description="Number of Spanning Reads">
##FORMAT=<ID=CGL,Number=R,Type=Integer,Description="Copies Gained or Lost with respect to reference">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$RUN_NAME
];

    print $spanN_vcffile $vcf_header;
    print $allwithsupport_vcffile $vcf_header;

    # Note: This query ~~simply~~ fetches all representative read TRs
    # for all supported alleles of all supported ref TRs, including
    # the reference allele (if supported).
    # This means ORDER is important. If changes to the query need
    # to be made, pay attention that the reference supporting read
    # is the FIRST one for each TR. Then the reads for the
    # remaining alleles, in increasing copy number.
    # Because the ordering of items in GROUP_CONCAT in SQLite is
    # not predictable, I use a subquery here instead to ensure
    # the ordering I want.
    my $get_supported_reftrs_query = qq{
        SELECT rid, is_singleton, firstindex, arlen, copynum,
          pattern AS consenuspat, head, REPLACE(UPPER(sequence), " ", "") AS sequence,
          refdir, GROUP_CONCAT(copies) AS cgl, (MAX(sameasref) == 1) AS refdetected,
          GROUP_CONCAT(support) AS support, GROUP_CONCAT(readarray) AS alt_seqs,
          GROUP_CONCAT(readdir) AS readdir, COUNT(*) AS num_alleles
        FROM (
            SELECT reftab.rid, is_singleton, firstindex, (lastindex - firstindex) + 1 AS arlen,
              reftab.copynum, reftab.pattern, reftab.head, sequence,
              c1.direction AS refdir, copies, sameasref, support,
              REPLACE(UPPER(SUBSTR(dna, first, (last-first)+1)), " ", "") AS readarray,
              c2.direction AS readdir
            FROM main.fasta_ref_reps mainreftab
              JOIN vntr_support ON mainreftab.rid = -vntr_support.refid
              JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
              JOIN clusterlnk c1 ON vntr_support.refid = c1.repeatid
              JOIN replnk ON vntr_support.representative = replnk.rid
              JOIN clusterlnk c2 ON c2.repeatid = replnk.rid
              JOIN fasta_reads ON replnk.sid = fasta_reads.sid
            WHERE support >= $MIN_SUPPORT_REQUIRED
            ORDER BY reftab.head ASC, reftab.firstindex ASC, sameasref DESC, copies ASC
        )
        GROUP BY rid};
    my $get_supported_reftrs_sth = $dbh->prepare($get_supported_reftrs_query);
    $get_supported_reftrs_sth->execute()
        or die "Cannot execute: " . $get_supported_reftrs_sth->errstr();

    # Bind columns from SELECT to hash
    my %row;
    $get_supported_reftrs_sth->bind_columns(
        \( @row{ @{ $get_supported_reftrs_sth->{NAME_lc} } } ) );

    # Loop over all supported ref TRs

    # If we're now seeing a new TR, write out the record for the
    # previous one (unless the rid is -1)
    my $vntr_count = 0;
    while ( $get_supported_reftrs_sth->fetch() ) {

        # If homozygous genotype called, either "0/0" or "1/1",
        # depending on whether or not the reference was detected.
        # Else, join a sequence from either 0 or 1 until the number
        # of alleles detected.
        my @gt
            = ( $row{num_alleles} == 1 )
            ? ( 1 * !$row{refdetected} ) x 2
            : ( 1 * !$row{refdetected}
                .. ( $row{num_alleles} - $row{refdetected} ) );

        # Split fields, and modify as needed
        my @alt_seqs     = split /,/, $row{alt_seqs};
        my @read_dirs    = split /,/, $row{readdir};
        my @cgl          = split /,/, $row{cgl};
        my @read_support = split /,/, $row{support};

        # If there is no DNA string for an allele, exit with error
        my @err_alleles = map { ( $alt_seqs[$_] eq "" ) ? ( $cgl[$_] ) : () }
            0 .. @cgl - 1;
        if (@err_alleles) {
            die "Error: sequence not found in database for ref ($row{rid}) alternate allele(s): ",
                join( ", ", @err_alleles ), "! Stopped at";
        }

        if ( $row{refdetected} ) {

            # Shift off ref sequence if reference was detected
            # (From note above, query result order matters)
            shift @alt_seqs;
            shift @read_dirs;
        }
        else {
            # New with 1.10, VCF 4.2 spec says FORMAT specs with
            # Number=R need to specify ref allele always, so
            # here we place a missing value (".") for it when we
            # don't see it for those fields
            unshift @cgl,          ".";
            unshift @read_support, ".";
        }

        # Reverse complement sequences if opposite dirs.
        # VCF spec says site with no alternate alleles gets
        # a missing value (".") for ALT field, so set default
        # to list containing just ".".
        my @rev_seqs = (".");
        @rev_seqs = map {
                  ( $read_dirs[$_] ne $row{refdir} )
                ? ( RC( $alt_seqs[$_] ) )
                : ( $alt_seqs[$_] )
            } 0 .. @read_dirs - 1
            if (@alt_seqs);

        $row{sequence}  = ( $row{sequence} ) ? $row{sequence} : ".";
        $row{support}   = join( ",", @read_support );
        $row{cgl}       = join( ",", @cgl );
        $row{alt}       = join( ",", @rev_seqs );
        $row{gt_string} = join( "/", @gt );

        my $vcf_rec = make_vcf_rec( \%row );

        # Only print record to spanN file if VNTR
        if ( $row{num_alleles} > 1 || !$row{refdetected} ) {
            $vntr_count++;
            print $spanN_vcffile $vcf_rec;
        }
        print $allwithsupport_vcffile $vcf_rec;
    }

    warn "Warning: Mismatch in VNTR count. Supported vntr count is $numvntrs"
        . " but we counted $vntr_count when producing VCF files."
        . " VCF header will have a bad VNTR count.\n"
        unless ( $vntr_count == $numvntrs );
    close $spanN_vcffile;
    close $allwithsupport_vcffile;
    $dbh->disconnect();
}

####################################

sub calc_entropy {
    my @ACGTcount = (0) x 4;
    my @diversity = ();
    my $in        = uc(shift);
    my $len       = length($in);
    my $i         = 0;
    my $count     = 0;

    while ( $i < $len ) {
        my $chr = substr( $in, $i, 1 );
        if ( $chr eq 'A' ) { $ACGTcount[0]++; $count++; }
        if ( $chr eq 'C' ) { $ACGTcount[1]++; $count++; }
        if ( $chr eq 'G' ) { $ACGTcount[2]++; $count++; }
        if ( $chr eq 'T' ) { $ACGTcount[3]++; $count++; }
        $i++;
    }

    if ( $count == 0 ) { return 0.0; }

    $diversity[0] = $ACGTcount[0] / $count;
    $diversity[1] = $ACGTcount[1] / $count;
    $diversity[2] = $ACGTcount[2] / $count;
    $diversity[3] = $ACGTcount[3] / $count;

    for ( my $e = 2; $e >= 0; $e-- ) {
        for ( my $f = 0; $f <= $e; $f++ ) {
            if ( $diversity[$f] < $diversity[ $f + 1 ] ) {
                my $temp = $diversity[ $f + 1 ];
                $diversity[ $f + 1 ] = $diversity[$f];
                $diversity[$f] = $temp;
            }
        }
    }

    my $entropy = (
        (     ( $diversity[0] == 0 ) ? 0
            : ( $diversity[0] * ( log( $diversity[0] ) / log(2) ) )
        ) + (
            ( $diversity[1] == 0 ) ? 0
            : ( $diversity[1] * ( log( $diversity[1] ) / log(2) ) )
        ) + (
            ( $diversity[2] == 0 ) ? 0
            : ( $diversity[2] * ( log( $diversity[2] ) / log(2) ) )
        ) + (
            ( $diversity[3] == 0 ) ? 0
            : ( $diversity[3] * ( log( $diversity[3] ) / log(2) ) )
        )
    );

    if ( $entropy < 0 ) { $entropy = -$entropy; }

    return $entropy;
}

############ Main #############

my $num;
my $i;

set_statistics( { N_MIN_SUPPORT => $MIN_SUPPORT_REQUIRED } );

my $dbh = get_dbh( { userefdb => 1 } )
    or die "Could not connect to database: $DBI::errstr";

#goto AAA;

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");

# Store temporary tables in memory
$dbh->do("PRAGMA temp_store = 2");

$dbh->do("DELETE FROM main.fasta_ref_reps")
    or die "Couldn't do statement: " . $dbh->errstr();

# Create a temporary table of all refids from vntr_support
# table, with the data aggregated to make parsing easier.
# Must negate refid to match in JOIN later.
$dbh->do(q{
    CREATE TEMPORARY TABLE temp.vntr_support AS
    SELECT rid, pattern, sameasref, support
    FROM refdb.fasta_ref_reps
      JOIN (
        SELECT -refid AS rid, GROUP_CONCAT(sameasref) AS sameasref,
          GROUP_CONCAT(support) AS support
        FROM vntr_support
        GROUP BY refid
        ORDER BY -refid ASC
    ) USING (rid)});

# Avoid DB lock issues with shared database
$dbh->do("DETACH DATABASE refdb");

my ($supported_vntr_count) = $dbh->selectrow_array(
    "SELECT COUNT(DISTINCT rid) FROM temp.vntr_support");

# now we can SELECT from reference TR table, using the temp table
# to JOIN.
my $get_supported_reftrs_sth = $dbh->prepare(
    "SELECT * FROM temp.vntr_support")
    or die "Couldn't prepare statement: " . $dbh->errstr();

print "Updating fasta_ref_reps table.\n";

my $query = q{
    INSERT INTO main.fasta_ref_reps (
        rid, alleles_sup, allele_sup_same_as_ref, entropy,
        has_support, span1, spanN, homez_same, homez_diff,
        hetez_same, hetez_diff, hetez_multi, support_vntr,
        support_vntr_span1 )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
};

my $update_ref_table_sth = $dbh->prepare($query)
    or die "Couldn't prepare statement: " . $dbh->errstr();

my $ReadTRsSupport = 0;
my @supported_refTRs;

# my $refid = -1;
$get_supported_reftrs_sth->execute()
    or die "Cannot execute: " . $get_supported_reftrs_sth->errstr();
$i = 0;
while ( my @data = $get_supported_reftrs_sth->fetchrow_array() ) {
    $i++;

    # Note: refid is positive
    my $ref = {
        refid              => $data[0],
        entr               => calc_entropy( $data[1] ),
        homez_same         => 0,
        homez_diff         => 0,
        hetez_same         => 0,
        hetez_diff         => 0,
        hetez_multi        => 0,
        support_vntr       => 0,
        span1              => 0,
        spanN              => 0,
        nsupport           => 0,
        has_support        => 0,
        nsameasref         => 0,
        readsum            => 0,
        support_vntr_span1 => 0
    };

    warn "TR $i; ", join( "; ", @data ), "\n"
        if ( $ENV{DEBUG} );
    my @sameasref = split ",", $data[2];
    my @support   = split ",", $data[3];

    for my $i ( 0 .. $#sameasref ) {
        if ( $support[$i] >= $MIN_SUPPORT_REQUIRED ) {

            # Count of alleles having at least MIN_SUPPORT_REQUIRED
            $ref->{nsupport}++;

            # Flag if ref has at least MIN_SUPPORT_REQUIRED
            ( $ref->{nsameasref} || $sameasref[$i] )
                && ( $ref->{nsameasref} = 1 );

            # Running sum of read TRs supporting alleles
            $ReadTRsSupport += $support[$i];
        }

        # Increments with support for each allele
        $ref->{readsum} += $support[$i];

        # Flag if this is a vntr with at least one read support
        ( $ref->{support_vntr_span1}
                || ( $support[$i] > 0 && $sameasref[$i] == 0 ) )
            && ( $ref->{support_vntr_span1} = 1 );
    }

    # Flag: at least one allele has at least MIN_SUPPORT_REQUIRED
    ( $ref->{nsupport} > 0 ) && ( $ref->{has_support} = 1 );

    ( $ref->{readsum} >= 1 )                     && ( $ref->{span1} = 1 );
    ( $ref->{readsum} >= $MIN_SUPPORT_REQUIRED ) && ( $ref->{spanN} = 1 );

    # hetez_multi is true if there is support
    # for more than one allele
    # Other genotype classes can only be true if hetez_multi
    # is false.
    # If supported alleles exceeds ploidy, call is multi.
    # Modify last branch to deal with special haploid case
    if ( $ref->{nsupport} > $run_conf{PLOIDY} ) {
        $ref->{hetez_multi} = 1;
    }
    elsif ( $ref->{nsupport} > 1 && $ref->{nsupport} <= $run_conf{PLOIDY} ) {
        ( $ref->{nsameasref} == 1 )
            ? ( $ref->{hetez_same} = 1 )
            : ( $ref->{hetez_diff} = 1 );
    }
    elsif ( $ref->{nsupport} == 1 ) {
        ( $ref->{nsameasref} == 1 )
            ? ( $ref->{homez_same} = 1 )
            : ( $ref->{homez_diff} = 1 );
    }

    $ref->{support_vntr}
        = 1 * ($ref->{homez_diff}
            || $ref->{hetez_same}
            || $ref->{hetez_diff}
            || $ref->{hetez_multi} );

    push(
        @supported_refTRs,

        # rid, alleles_sup, allele_sup_same_as_ref, entropy,
        #     has_support, span1, spanN, homez_same, homez_diff,
        #     hetez_same, hetez_diff, hetez_multi, support_vntr,
        #     support_vntr_span1
        [   $ref->{refid},        $ref->{nsupport},
            $ref->{nsameasref},   $ref->{entr},
            $ref->{has_support},  $ref->{span1},
            $ref->{spanN},        $ref->{homez_same},
            $ref->{homez_diff},   $ref->{hetez_same},
            $ref->{hetez_diff},   $ref->{hetez_multi},
            $ref->{support_vntr}, $ref->{support_vntr_span1}
        ]
    );

    if ( @supported_refTRs % $RECORDS_PER_INFILE_INSERT == 0 ) {
        my $cb   = gen_exec_array_cb( \@supported_refTRs );
        my $rows = vs_db_insert( $dbh, $update_ref_table_sth, $cb,
            "Error when inserting entries into our supported reference TRs table.");
        @supported_refTRs = ();
    }
}

my $updfromtable = $i;

# Insert last rows:
if (@supported_refTRs) {
    my $cb   = gen_exec_array_cb( \@supported_refTRs );
    my $rows = vs_db_insert( $dbh, $update_ref_table_sth, $cb,
        "Error when inserting entries into our supported reference TRs table.");
    @supported_refTRs = ();
}
if ( $updfromtable != $supported_vntr_count ) {
    die "Updated number of entries($updfromtable) not equal to the number of references, aborting!";
}

# updating stats table
print "Updating stats table.\n";
my $update_stats_sth = $dbh->prepare(q{
    UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED = (
        SELECT COUNT(*)
        FROM (
            SELECT COUNT(*) AS thecount
            FROM clusterlnk
            WHERE repeatid < 0
            GROUP BY clusterid HAVING thecount = 1
        ) f
    )})
    or die "Couldn't prepare statement: " . $dbh->errstr();
$update_stats_sth->execute()
    or die "Cannot execute: " . $update_stats_sth->errstr();

$dbh->do('CREATE TEMPORARY TABLE t1 (c1 INT PRIMARY KEY NOT NULL)')
    or die "Couldn't do statement: " . $dbh->errstr();

my $sth = $dbh->prepare(q{
    INSERT INTO t1
    SELECT urefid
    FROM (
        SELECT COUNT(*) AS thecount, MAX(repeatid) AS urefid
        FROM clusterlnk
        WHERE repeatid < 0
        GROUP BY clusterid
        HAVING thecount = 1
    ) f})
    or die "Couldn't prepare statement: " . $dbh->errstr();
$sth->execute() or die "Cannot execute: " . $sth->errstr();

$sth = $dbh->prepare(q{
    UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER = (
        SELECT COUNT(DISTINCT repeatid)
        FROM clusterlnk
        WHERE repeatid IN (SELECT c1 FROM t1)
    )})
    or die "Couldn't prepare statement: " . $dbh->errstr();
$sth->execute() or die "Cannot execute: " . $sth->errstr();

$sth = $dbh->prepare(q{
    UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED = (
        SELECT COUNT(distinct map.refid)
        FROM clusterlnk
          JOIN map ON map.refid = -clusterlnk.repeatid
        WHERE repeatid IN (select c1 from t1)
    )})
    or die "Couldn't prepare statement: " . $dbh->errstr();
$sth->execute() or die "Cannot execute: " . $sth->errstr();


$sth = $dbh->prepare(q{
    UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED =
        NUMBER_REFS_SINGLE_REF_CLUSTER -
        NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED})
    or die "Couldn't prepare statement: " . $dbh->errstr();
$sth->execute() or die "Cannot execute: " . $sth->errstr();

# Update last few stats:
my ( $mapped, $rank, $rankflank, $num_spanN );

($mapped) = $dbh->selectrow_array("SELECT COUNT(*) FROM map")
    or die "Couldn't select map count: " . $dbh->errstr();

($rank) = $dbh->selectrow_array("SELECT COUNT(*) FROM rank")
    or die "Couldn't select rank count: " . $dbh->errstr();

($rankflank) = $dbh->selectrow_array("SELECT COUNT(*) FROM rankflank")
    or die "Couldn't select rankflank count: " . $dbh->errstr();

# update spanN number on stats
($num_spanN) = $dbh->selectrow_array(q{
    SELECT COUNT(*)
    FROM fasta_ref_reps
    WHERE support_vntr > 0})
    or die "Couldn't select span N count: " . $dbh->errstr();

$dbh->do(qq{
    UPDATE stats SET
      NUMBER_MAPPED = $mapped,
      NUMBER_RANK = $rank,
      NUMBER_RANKFLANK = $rankflank,
      NUMBER_REFS_VNTR_SPAN_N = $num_spanN})
      or die "Couldn't do statement: " . $dbh->errstr();

# set old db settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");
$dbh->do("PRAGMA temp_store = 0");

# Optimize queries
$dbh->do("PRAGMA main.optimize");
$dbh->disconnect();

# Print VCF
print "Producing VCFs.\n";
print_vcf();
