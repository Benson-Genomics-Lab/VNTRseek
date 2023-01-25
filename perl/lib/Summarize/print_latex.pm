package Summarize::print_latex;

# Kyler Anderson 1/5/23
# Code from a previous author that's not currently executed
# Original use in updaterefs.pl

use v5.24;
use carp;

use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil qw(get_dbh get_statistics);

#use GD::Graph::linespoints;

use base 'Exporter';
our @EXPORT_OK = 'print_latex';

####################################

sub commify {
    my $input = shift;
    warn "Undef input" unless defined $input;
    $input = reverse $input;
    $input =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return reverse $input;
}

###################

sub print_latex {
    if (@_ < 3) { carp "print_dist expects 3 args"; }
    my $ReadTRsSupport = $_[0];
    my $result_prefix = $_[1];
    my $MIN_SUPPORT_REQUIRED = $_[2];

    # warn "Read TRs supported: $ReadTRsSupport\n"
    #     if ($ENV{DEBUG});

    my $sum_has_support  = 0;
    my $sum_span1        = 0;
    my $sum_spanN        = 0;
    my $sum_homez_same   = 0;
    my $sum_homez_diff   = 0;
    my $sum_hetez_same   = 0;
    my $sum_hetez_diff   = 0;
    my $sum_hetez_multi  = 0;
    my $sum_support_vntr = 0;

    # Get needed stats
    my @stats = qw(NUMBER_READS
        NUMBER_TRS_IN_READS
        NUMBER_TRS_IN_READS_GE7
        NUMBER_READS_WITHTRS_GE7
        NUMBER_READS_WITHTRS
        NUMBER_READS_WITHTRS_GE7_AFTER_REDUND
        CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS
        NUMBER_REF_TRS
        NUMBER_REFS_TRS_AFTER_REDUND
        CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS
        CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR);
    my $stat_hash = get_statistics(@stats);

    my $dbh = get_dbh( { readonly => 1, userefdb => 1 } );
    my $sth = $dbh->prepare(
        q{SELECT sum(has_support),sum(span1),sum(spanN),sum(homez_same),sum(homez_diff),
            sum(hetez_same),sum(hetez_diff),sum(hetez_multi),sum(support_vntr)
            FROM main.fasta_ref_reps}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {

        $sum_has_support  = $data[0];
        $sum_span1        = $data[1];
        $sum_spanN        = $data[2];
        $sum_homez_same   = $data[3];
        $sum_homez_diff   = $data[4];
        $sum_hetez_same   = $data[5];
        $sum_hetez_diff   = $data[6];
        $sum_hetez_multi  = $data[7];
        $sum_support_vntr = $data[8];
    }
    $sth->finish();

# Definitions
# READS -- Reads from the subject data set
# READ-TRs -- read TRs from the subject data set
# REF-TRs -- reference TRs from the hg19 data set
# INITIAL -- total started with in the original data set
# GE7 -- only patterns >= 7
# RDE (Redundancy Elimination) -- left after redundancy elimination
# CRDE (Cyclic Redundancy Elimination) -- left after cyclic redundancy elimination
# PC (Profile Clustered) -- PROCLU clustered into a cluster with at least one reference and at least one read
# MAP (Mapped) -- flank TTT clustered into a cluster with one reference and at least one read
# TIE-OK -- applies to a read which is MAP, but may appear in more than one flank TTT cluster
#    (due to ties or no flanking sequence)
# ADDBACK -- includes those eliminated by CRDE, but added back for TTT flank alignment
# SINGLETON -- PC into a cluster with EXACTLY ONE REFERENCE and MAP
# DISTINGUISHABLE -- PC into a cluster with MORE THAN ONE REFERENCE but flank distinguishable and MAP
# INDISTINGUISHABLE -- PC into a cluster with MORE THAN ONE REFERENCE, NOT flank indistinguishable and MAP
# SPAN1 -- MAP with at least 1 read that is NOT assembly required
# SPANN -- MAP with at least N reads that are NOT assembly required
# SUPPORT -- MAP with support equal to 2 or more for at least one copy number
# HOMOZYGOUS -- MAP only one copy number has SUPPORT
# HETEROZYGOUS -- MAP two or more copy numbers have SUPPORT
# MULTI -- MAP three or more copy numbers have support
# SAME -- one copy number with SUPPORT matches reference copy number
# DIFF -- NO copy number with SUPPORT matches reference copy number

    # ***************************************************
    # Yevgeniy, assign real numbers to these variables.

    my $totalReads = $stat_hash->{NUMBER_READS};    # INITIAL READS

#$totalReadsWithTRs = $data[0];       # INITIAL READS which contain TRs
#$totalReadsWithTRsPatternGE7  = $data[0];        # INITIAL READS which contain TRs GE7
#$readTRsWithPatternGE7  = $data[0];       # INITIAL READ-TRs GE7
#$readTRsWPGE7AfterCyclicRedundancyElimination = $data[0];       # INITIAL READ-TRs GE7 CRDE

    my $totalReadTRs = $stat_hash->{NUMBER_TRS_IN_READS};

    my $readTRsWithPatternGE7       = $stat_hash->{NUMBER_TRS_IN_READS_GE7};
    my $totalReadsWithTRsPatternGE7 = $stat_hash->{NUMBER_READS_WITHTRS_GE7};
    my $totalReadsWithTRs           = $stat_hash->{NUMBER_READS_WITHTRS};
    my $readTRsWPGE7AfterCyclicRedundancyElimination
        = $stat_hash->{NUMBER_READS_WITHTRS_GE7_AFTER_REDUND};

    my $readTRsProfileClustered
        = $stat_hash->{CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS}; 
        # INITIAL READ-TRs GE7 PC ADDBACK

    $sth
        = $dbh->prepare(
        "SELECT count(distinct refid), count(distinct readid) FROM map WHERE bbb=1"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    my $refTRsMapped  = 0;
    my $readTRsMapped = 0;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    {
        my @data = $sth->fetchrow_array();
        $refTRsMapped = $data[0];

        # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP TIE-OK
        $readTRsMapped = $data[1];
    }
    $sth->finish();

    my $refTRsAREWithPatternGE7
        = $stat_hash->{NUMBER_REF_TRS};    # INITIAL REF-TRs RDE GE7
    my $refTRsAREWPGE7AfterCyclicRedundancyElimination
        = $stat_hash->{NUMBER_REFS_TRS_AFTER_REDUND}
        ;                                  # INITIAL REF-TRs RDE GE7 CRDE
    my $refTRsProfileClustered
        = $stat_hash->{CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS}
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK

    my $refTRsMappedSpan1
        = $sum_span1;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SPAN1
    my $refTRsMappedSpanN
        = $sum_spanN;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SPANN
    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SINGLETON
    #my $refTRsMappedSingleton = GetStatistics("NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED");

    # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP SINGLETON
    my $readTRsMappedToSingleton = 0;

    # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP TIE-OK INDISTINGUISHABLE
    my $readTRsMappedToIndistinguishable = 0;

    $sth = $dbh->prepare(
        q{SELECT is_singleton, count(distinct map.readid)
    FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid
    WHERE bbb=1 GROUP BY is_singleton}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        ( $data[0] == 1 ) && ( $readTRsMappedToSingleton         = $data[1] );
        ( $data[0] == 0 ) && ( $readTRsMappedToIndistinguishable = $data[1] );
    }

    unless ( $readTRsMapped
        == $readTRsMappedToIndistinguishable + $readTRsMappedToSingleton ) {
        die "Error: mismatch of sum of mapped read TRs by distinguishablility and total read TRs mapped. "
            . "(Expected $readTRsMapped, got"
            . " $readTRsMappedToIndistinguishable + $readTRsMappedToSingleton = "
            . ($readTRsMappedToIndistinguishable + $readTRsMappedToSingleton )
            . ")\n";
    }

    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SINGLETON
    my $refTRsMappedSingleton = 0;

    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP INDISTINGUISHABLE
    my $refTRsMappedIndistinguishable   = 0;
    my $InvarTRsMappedSingleton         = 0;
    my $InvarTRsMappedIndistinguishable = 0;

    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SINGLETON VNTR
    my $VNTRasSingleton = 0;

    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP INDISTINGUISHABLE VNTR
    my $VNTRasIndistinguishable = 0;

    $sth = $dbh->prepare(
        q{SELECT is_singleton, support_vntr, count(distinct rid)
    FROM main.fasta_ref_reps mainref INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
    INNER JOIN map ON mainref.rid=map.refid
    WHERE bbb=1 GROUP BY is_singleton, support_vntr}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        ( $data[0] == 0 )
            && ( $data[1] == 0 )
            && ( $InvarTRsMappedIndistinguishable = $data[2] );
        ( $data[0] == 0 )
            && ( $data[1] == 1 )
            && ( $VNTRasIndistinguishable = $data[2] );
        ( $data[0] == 1 )
            && ( $data[1] == 0 )
            && ( $InvarTRsMappedSingleton = $data[2] );
        ( $data[0] == 1 )
            && ( $data[1] == 1 )
            && ( $VNTRasSingleton = $data[2] );
    }

    $refTRsMappedSingleton = $InvarTRsMappedSingleton + $VNTRasSingleton;
    $refTRsMappedIndistinguishable
        = $InvarTRsMappedIndistinguishable + $VNTRasIndistinguishable;

    unless ( $refTRsMapped == $refTRsMappedIndistinguishable + $refTRsMappedSingleton ) {
        die "Error: mismatch of sum of mapped ref TRs by distinguishablility and total ref TRs mapped. "
            . "(Expected $refTRsMapped, got '$refTRsMappedIndistinguishable' + '$refTRsMappedSingleton' = '"
            . ( $refTRsMappedIndistinguishable + $refTRsMappedSingleton )
            . "')\n";
    }

    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT
    my $refWithSupport = $sum_has_support;
    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HOMZYGOUS SAME
    my $refAsHomozygousSame = $sum_homez_same;
    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HOMZYGOUS DIFF
    my $refAsHomozygousDiff = $sum_homez_diff;
    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HETEROZYGOUS SAME
    my $refAsHeterozygousSame = $sum_hetez_same;
    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HETEROZYGOUS DIFF
    my $refAsHeterozygousDiff = $sum_hetez_diff;
    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT MULTI
    my $refAsMulti = $sum_hetez_multi;
    # ***************************************************

    # formatted to show commas at thousands
    my $form_totalReads        = commify($totalReads);
    my $form_totalReadsWithTRs = commify($totalReadsWithTRs);
    my $form_totalReadsWithTRsPatternGE7
        = commify($totalReadsWithTRsPatternGE7);
    my $form_totalReadTRs          = commify($totalReadTRs);
    my $form_readTRsWithPatternGE7 = commify($readTRsWithPatternGE7);
    my $form_readTRsWPGE7AfterCyclicRedundancyElimination
        = commify($readTRsWPGE7AfterCyclicRedundancyElimination);
    my $form_readTRsProfileClustered  = commify($readTRsProfileClustered);
    my $form_readTRsMapped            = commify($readTRsMapped);
    my $form_readTRsMappedToSingleton = commify($readTRsMappedToSingleton);
    my $form_readTRsMappedToIndistinguishable
        = commify($readTRsMappedToIndistinguishable);

    my $form_refTRsAREWithPatternGE7 = commify($refTRsAREWithPatternGE7);
    my $form_refTRsAREWPGE7AfterCyclicRedundancyElimination
        = commify($refTRsAREWPGE7AfterCyclicRedundancyElimination);
    my $form_refTRsProfileClustered = commify($refTRsProfileClustered);
    my $form_refTRsMapped           = commify($refTRsMapped);
    my $form_refTRsMappedSpan1      = commify($refTRsMappedSpan1);
    my $form_refTRsMappedSpanN      = commify($refTRsMappedSpanN);
    my $form_refTRsMappedSingleton  = commify($refTRsMappedSingleton);
    my $form_refTRsMappedIndistinguishable
        = commify($refTRsMappedIndistinguishable);

    my $form_VNTRasSingleton         = commify($VNTRasSingleton);
    my $form_VNTRasIndistinguishable = commify($VNTRasIndistinguishable);

    my $form_refWithSupport        = commify($refWithSupport);
    my $form_refAsHomozygousSame   = commify($refAsHomozygousSame);
    my $form_refAsHomozygousDiff   = commify($refAsHomozygousDiff);
    my $form_refAsHeterozygousSame = commify($refAsHeterozygousSame);
    my $form_refAsHeterozygousDiff = commify($refAsHeterozygousDiff);
    my $form_refAsMulti            = commify($refAsMulti);
    my $form_ReadTRsSupport        = commify($ReadTRsSupport);

    # DO NOT COMPUTE
    my $percentReadTRsProfileClustered;
    my $percentReadTRsMapped;
    my $percentRefTRsProfileClustered;
    my $percentRefTRsMapped;
    my $percentRefTRsMappedSingleton;
    my $percentRefTRsMappedIndistinguishable;
    my $percentReadTRsMappedToSingleton;
    my $percentReadTRsMappedToIndistinguishable;
    my $percentVNTRasSingleton;
    my $percentVNTRasIndistinguishable;
    my $percentRefTRsMappedSpan1;
    my $percentRefTRsMappedSpanN;
    my $percentRefWithSupport;
    my $percentRefAsHomozygousSame;
    my $percentRefAsHomozygousDiff;
    my $percentRefAsHeterozygousSame;
    my $percentRefAsHeterozygousDiff;
    my $percentRefAsMulti;
    my $percentVNTRTotalByGeneticClass;
    my $VNTRTotalByRefClass
        = commify( $stat_hash->{CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR} );
    my $VNTRTotalByGeneticClass;

    # latex header
    open( my $texfh, ">", "${result_prefix}.span${MIN_SUPPORT_REQUIRED}.tex" )
        or die
        "\nCan't open for writing ${result_prefix}.span${MIN_SUPPORT_REQUIRED}.tex\n\n";
    printf( $texfh "\n\\documentclass[twoside,10pt]{article}"
            . "\n\\usepackage{underscore}"
            . "\n\\setlength{\\topmargin}{-1cm}"
            . "\n\\setlength{\\oddsidemargin}{-0.5cm}"
            . "\n\\setlength{\\evensidemargin}{-0.5cm}"
            . "\n\\setlength{\\textwidth}{6.5in}"
            . "\n\\setlength{\\textheight}{9in}"
            . "\n\\setlength{\\parskip}{0.1 in}"
            . "\n\\setlength{\\parindent}{0.0 in}"
            . "\n\\begin{document}"

            . "\n\\begin{center}"
            . "\n{\\Large{\\bf VNTRseek analysis of %s}}"
            . "\n\\end{center}", $RUN_NAME
    );

    # TR Clustering and Mapping Outcomes
    $percentRefTRsProfileClustered
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsProfileClustered / $refTRsAREWithPatternGE7 )
        : 0;
    $percentRefTRsMapped
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsMapped / $refTRsAREWithPatternGE7 )
        : 0;
    $percentReadTRsProfileClustered
        = $readTRsWithPatternGE7
        ? int( 100 * $readTRsProfileClustered / $readTRsWithPatternGE7 )
        : 0;
    $percentReadTRsMapped
        = $readTRsWithPatternGE7
        ? int( 100 * $readTRsMapped / $readTRsWithPatternGE7 )
        : 0;
    $percentRefTRsMappedSingleton
        = $refTRsMapped
        ? int( 100 * $refTRsMappedSingleton / $refTRsMapped )
        : 0;
    $percentRefTRsMappedIndistinguishable
        = $refTRsMapped
        ? int( 100 * $refTRsMappedIndistinguishable / $refTRsMapped )
        : 0;
    $percentReadTRsMappedToSingleton
        = $readTRsMapped
        ? int( 100 * $readTRsMappedToSingleton / $readTRsMapped )
        : 0;
    $percentReadTRsMappedToIndistinguishable
        = $readTRsMapped
        ? int( 100 * $readTRsMappedToIndistinguishable / $readTRsMapped )
        : 0;
    $VNTRTotalByRefClass = $VNTRasSingleton + $VNTRasIndistinguishable;
    my $form_VNTRTotalByRefClass = commify($VNTRTotalByRefClass);
    $percentVNTRasSingleton
        = $VNTRTotalByRefClass
        ? int( 100 * $VNTRasSingleton / $VNTRTotalByRefClass )
        : 0;
    $percentVNTRasIndistinguishable
        = $VNTRTotalByRefClass
        ? int( 100 * $VNTRasIndistinguishable / $VNTRTotalByRefClass )
        : 0;

    # VNTR calling
    $percentRefTRsMappedSpan1
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsMappedSpan1 / $refTRsAREWithPatternGE7 )
        : 0;
    $percentRefTRsMappedSpanN
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsMappedSpanN / $refTRsAREWithPatternGE7 )
        : 0;
    $percentRefWithSupport
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refWithSupport / $refTRsAREWithPatternGE7 )
        : 0;
    $VNTRTotalByGeneticClass
        = $refAsHomozygousDiff
        + $refAsHeterozygousSame
        + $refAsHeterozygousDiff
        + $refAsMulti;

    #$form_VNTRTotalByGeneticClass = commify($VNTRTotalByGeneticClass);
    $percentRefAsHomozygousSame
        = $refWithSupport
        ? int( 100 * $refAsHomozygousSame / $refWithSupport )
        : 0;
    $percentRefAsHomozygousDiff
        = $refWithSupport
        ? int( 100 * $refAsHomozygousDiff / $refWithSupport )
        : 0;
    $percentRefAsHeterozygousSame
        = $refWithSupport
        ? int( 100 * $refAsHeterozygousSame / $refWithSupport )
        : 0;
    $percentRefAsHeterozygousDiff
        = $refWithSupport
        ? int( 100 * $refAsHeterozygousDiff / $refWithSupport )
        : 0;
    $percentRefAsMulti
        = $refWithSupport ? int( 100 * $refAsMulti / $refWithSupport ) : 0;
    $percentVNTRTotalByGeneticClass
        = $refWithSupport
        ? int( 100 * $VNTRTotalByGeneticClass / $refWithSupport )
        : 0;

#$percentReadTRsSupport = $readTRsWithPatternGE7 ? int(100*$ReadTRsSupport/$readTRsWithPatternGE7) : 0;

    my $form_readsMapped   = "";
    my $percentReadsMapped = $form_readsMapped;

    $sth = $dbh->prepare(q{
            SELECT count(distinct head)
            FROM map
              JOIN replnk ON map.readid = replnk.rid
              JOIN fasta_reads ON fasta_reads.sid = replnk.sid
            WHERE bbb = 1})
            or die "Couldn't prepare statement: " . $dbh->errstr();
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    {
        my @data = $sth->fetchrow_array();
        $form_readsMapped = commify( int( $data[0] ) );
        $percentReadsMapped
            = $totalReadsWithTRsPatternGE7
            ? int( 100 * int( $data[0] ) / $totalReadsWithTRsPatternGE7 )
            : 0;
    }
    $sth->finish();

    #New Table Formats
    my $percentRefTRsMappedSpan1byMapped;
    my $percentRefTRsMappedSpanNbyMapped;
    my $percentRefWithSupportbyMapped;
    my $form_refOneAlleleDiff;
    my $form_refTwoAllelesSame;
    my $form_refTwoAllelesDiff;
    my $form_refMultiAlleles;
    my $percentRefOneAlleleDiff;
    my $percentRefTwoAllelesSame;
    my $percentRefTwoAllelesDiff;
    my $percentRefMultiAlleles;

    $percentRefTRsMappedSpan1byMapped
        = $refTRsMapped ? int( 100 * $refTRsMappedSpan1 / $refTRsMapped ) : 0;
    $percentRefTRsMappedSpanNbyMapped
        = $refTRsMapped ? int( 100 * $refTRsMappedSpanN / $refTRsMapped ) : 0;
    $percentRefWithSupportbyMapped
        = $refTRsMapped ? int( 100 * $refWithSupport / $refTRsMapped ) : 0;
    $form_refOneAlleleDiff  = commify($refAsHomozygousDiff);
    $form_refTwoAllelesSame = commify($refAsHeterozygousSame);
    $form_refTwoAllelesDiff = commify($refAsHeterozygousDiff);
    $form_refMultiAlleles   = commify($refAsMulti);
    $percentRefOneAlleleDiff
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsHomozygousDiff / $VNTRTotalByGeneticClass )
        : 0;
    $percentRefTwoAllelesSame
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsHeterozygousSame / $VNTRTotalByGeneticClass )
        : 0;
    $percentRefTwoAllelesDiff
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsHeterozygousDiff / $VNTRTotalByGeneticClass )
        : 0;
    $percentRefMultiAlleles
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsMulti / $VNTRTotalByGeneticClass )
        : 0;

    printf( $texfh "\n\\begin{table}[htdp]"
            . "\n\\begin{center}"
            . "\nA. Mapping\\\\"
            . "\n\\vspace{.1in}"
            . "\n\\begin{tabular}{|c|c||c|c|c||}"
            . "\n\\hline"
            . "\n\\multicolumn{5}{|c||}{Mapping}\\\\"
            . "\n\\hline"
            . "\n&&Input&&\\\\"
            . "\n&Total&(After Filters)&Mapped&\\%%\\\\"
            . "\n\\hline"
            . "\nReference TRs&---&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\nRead TRs&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\nReads&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\n\\end{tabular}\\\\"
            . "\n\\vspace{.2in}" . "\n"
            . "\nB. Reference Results\\\\"
            . "\n\\vspace{.1in}"
            . "\n\\begin{tabular}{|c|c||c||c|c||}"
            . "\n\\hline"
            . "\n\\multicolumn{5}{|c||}{Mapped Reference TRs}\\\\"
            . "\n\\hline"
            . "\n\\multicolumn{2}{|c||}{}&At Least&\\multicolumn{2}{c||}{}\\\\"
            . "\n\\multicolumn{2}{|c||}{Mapped by at least}&One Allele&\\multicolumn{2}{c||}{By Reference Category}\\\\"
            . "\n\\cline{1-2}"
            . "\n\\cline{4-5}"
            . "\nOne Read&Two Reads& Supported&Singleton&Indistinguishable\\\\\\hline"
            . "\n%s&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\n%s\\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%\\\\"
            . "\n\\hline"
            . "\n\\end{tabular}"
            . "\n\\vspace{.2in}" . "\n" . "\n"
            . "\nC. VNTR Results\\\\"
            . "\n\\vspace{.1in}"
            . "\n\\begin{tabular}{|c||c|c|c|c||c|c||c||}"
            . "\n\\hline"
            . "\n\\multicolumn{7}{|c||}{VNTRs}\\\\"
            . "\n\\hline"
            . "\n&\\multicolumn{4}{c||}{Alleles Supported}&\\multicolumn{2}{c||}{}\\\\"
            . "\n\\cline{2-5}"
            . "\n&One&\\multicolumn{3}{c||}{Two or More}&\\multicolumn{2}{c||}{By Reference Category}\\\\"
            . "\n\\cline{2-5}"
            . "\n&\$\\star\$&\$\\bullet\$&\$\\bullet\$&\$\\bullet\$&\\multicolumn{2}{c||}{}\\\\"
            . "\n\\cline{6-7}"
            . "\nTotal&Diff&Same/Diff&Diff/Diff &Multi&Singleton&Indistinguishable\\\\"
            . "\n\\hline"
            . "\n%s&%s&%s&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\n100\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%\\\\"
            . "\n\\hline"
            . "\n\\end{tabular}\\\\"
            . "\n\\ \\\\"
            . "\n\$\\star\$ Inferred VNTR \\quad \$\\bullet\$ Observed VNTR"
            . "\n\\vspace{.2in}"
            . "\n\\end{center}"
            . "\n\\caption{{\\bf VNTRseek Results.} }"
            . "\n\\label{table:mapping and mapped by reference category}"
            . "\n\\end{table}\%%" . "\n",
        $form_refTRsAREWithPatternGE7,
        $form_refTRsMapped,
        $percentRefTRsMapped,
        $form_totalReadTRs,
        $form_readTRsWithPatternGE7,
        $form_readTRsMapped,
        $percentReadTRsMapped,
        $form_totalReads,
        $form_totalReadsWithTRsPatternGE7,
        $form_readsMapped,
        $percentReadsMapped,
        $form_refTRsMappedSpan1,
        $form_refTRsMappedSpanN,
        $form_refWithSupport,
        $form_refTRsMappedSingleton,
        $form_refTRsMappedIndistinguishable,
        ,
        $percentRefTRsMappedSpan1byMapped,
        $percentRefTRsMappedSpanNbyMapped,
        $percentRefWithSupportbyMapped,
        $percentRefTRsMappedSingleton,
        $percentRefTRsMappedIndistinguishable,
        $form_VNTRTotalByRefClass,
        $form_refOneAlleleDiff,
        $form_refTwoAllelesSame,
        $form_refTwoAllelesDiff,
        $form_refMultiAlleles,
        $form_VNTRasSingleton,
        $form_VNTRasIndistinguishable,
        ,
        $percentRefOneAlleleDiff,
        $percentRefTwoAllelesSame,
        $percentRefTwoAllelesDiff,
        $percentRefMultiAlleles,
        $percentVNTRasSingleton,
        $percentVNTRasIndistinguishable
    );
    printf( $texfh "\n\\end{document}" );
    $dbh->disconnect;
    close($texfh);
}