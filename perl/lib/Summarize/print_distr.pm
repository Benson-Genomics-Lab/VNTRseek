package Summarize::print_distr;

# Kyler Anderson 1/5/23
# Code from a previous author that's not currently executed
# Original use in updaterefs.pl

use v5.24;
use carp;

use FindBin;
use lib "$FindBin::RealBin/lib";
use vutil 'get_dbh';

#use GD::Graph::linespoints;

use base 'Exporter';
our @EXPORT_OK = 'print_distr';

####################################

sub commify {
    my $input = shift;
    warn "Undef input" unless defined $input;
    $input = reverse $input;
    $input =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return reverse $input;
}

####################################

sub print_distr {
    if (@_ < 2) { carp "print_dist expects 2 args"; }
    my $result_prefix = $_[0];
    my $MIN_SUPPORT_REQUIRED = $_[1];

    my $LARGEST_PSIZE = 122;
    my $LARGEST_ASIZE = 311;

    #my $LARGEST_PSIZE = 5000;
    #my $LARGEST_ASIZE = 5000;

    my $BUCKET_PAT = 1;
    my $BUCKET_AS  = 20;

    my %PHASH1      = ();
    my %PHASH2      = ();
    my %PHASH3      = ();
    my %PHASH4      = ();
    my %PHASH5      = ();
    my %PHASH6      = ();
    my %PHASHB      = ();
    my %PHASH_SN    = ();
    my %PHASH_SN_VN = ();

    open my $distrfh, ">", "${result_prefix}.span${MIN_SUPPORT_REQUIRED}.txt"
        or die "\nCan't open for reading ${result_prefix}.span${MIN_SUPPORT_REQUIRED}.txt\n";
    print $distrfh "\n\nPatternSize, All, Span1, PercentageS1,"
        . " Span${MIN_SUPPORT_REQUIRED}, PercentageS${MIN_SUPPORT_REQUIRED},"
        . " VntrSpan${MIN_SUPPORT_REQUIRED}, PercentageVS${MIN_SUPPORT_REQUIRED}\n";

    my $dbh = get_dbh( { readonly => 1, userefdb => 1 } );

    # 1 patsize (all)
    my $sth = $dbh->prepare(q{
        SELECT length(pattern) AS patsize, count(*)
        FROM refdb.fasta_ref_reps
        GROUP BY patsize
        ORDER BY patsize ASC})
        or die "Couldn't prepare statement: " . $dbh->errstr();
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    my $i = 0;
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH1{ $data[0] } = $data[1];
        }
        else {
            $PHASH1{$LARGEST_PSIZE} += $data[1];
        }
        $i++;
    }

    # 2 patsize (span1)
    $sth = $dbh->prepare(q{
        SELECT length(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab
          JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.span1 > 0
        GROUP BY patsize
        ORDER BY patsize ASC})
        or die "Couldn't prepare statement: " . $dbh->errstr();
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH2{ $data[0] } = $data[1];
        }
        else {
            $PHASH2{$LARGEST_PSIZE} += $data[1];
        }
        $i++;
    }

    # 3 patsize (spanN)
    $sth = $dbh->prepare(q{
        SELECT length(pattern) AS patsize, count(*)
        FROM main.fasta_ref_reps mainreftab
          JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.spanN > 0
        GROUP BY patsize
        ORDER BY patsize ASC})
        or die "Couldn't prepare statement: " . $dbh->errstr();
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH_SN{ $data[0] } = $data[1];
        }
        else {
            $PHASH_SN{$LARGEST_PSIZE} += $data[1];
        }
        $i++;
    }

    # 4 patsize (vntr spanN)
    $sth = $dbh->prepare(q{
        SELECT LENGTH(pattern) AS patsize, count(*)
        FROM main.fasta_ref_reps mainreftab
          JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE support_vntr > 0
        GROUP BY patsize
        ORDER BY patsize ASC})
        or die "Couldn't prepare statement: " . $dbh->errstr();
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH_SN_VN{ $data[0] } = $data[1];
        }
        else {
            $PHASH_SN_VN{$LARGEST_PSIZE} += $data[1];
        }
        $i++;
    }

    # print
    my $maxval     = 10;
    my @arX0       = ();
    my @arX        = ();
    my @arY1       = ();
    my @arY2       = ();
    my %KEYS12     = ();
    my $tall       = 0;
    my $tspan1     = 0;
    my $tspanN     = 0;
    my $tspanNvntr = 0;

    foreach my $key ( keys %PHASH1 )      { $KEYS12{$key} = 1; }
    foreach my $key ( keys %PHASH2 )      { $KEYS12{$key} = 1; }
    foreach my $key ( keys %PHASH_SN )    { $KEYS12{$key} = 1; }
    foreach my $key ( keys %PHASH_SN_VN ) { $KEYS12{$key} = 1; }
    foreach my $key ( sort { $a <=> $b } ( keys %KEYS12 ) ) {
        my $val1 = 0;
        if ( exists $PHASH1{$key} ) { $val1 = $PHASH1{$key}; }
        my $val2 = 0;
        if ( exists $PHASH2{$key} ) { $val2 = $PHASH2{$key}; }
        my $val3 = 0;
        if ( exists $PHASH_SN{$key} ) { $val3 = $PHASH_SN{$key}; }
        my $val4 = 0;
        if ( exists $PHASH_SN_VN{$key} ) { $val4 = $PHASH_SN_VN{$key}; }

        if ( $key < $LARGEST_PSIZE ) {
            print $distrfh $key . ", "
                . $val1 . ", "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "\%\n";
            push( @arX0, "$key" );
        }
        else {
            print $distrfh $key . "+, "
                . $val1 . ", "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "\%\n";
            push( @arX0, "$key+" );
        }
        push( @arX,  $key );
        push( @arY1, $val1 );
        push( @arY2, $val2 );

        $tall       += $val1;
        $tspan1     += $val2;
        $tspanN     += $val3;
        $tspanNvntr += $val4;

        $maxval = max( $maxval, $val1 );
    }
    print "TOTAL, $tall, $tspan1, \n";
    print $distrfh "TOTAL, $tall, $tspan1, "
        . int( 100 * $tspan1 / $tall )
        . "\%, $tspanN, "
        . int( 100 * $tspanN / $tall )
        . "\%, $tspanNvntr, "
        . int( 100 * $tspanNvntr / $tall ) . "\%\n";

#my $temp = commify(GetStatistics("NUMBER_REF_TRS"));
#(my $day, my $month, my $year) = (localtime)[3,4,5];
#my $temp2 = sprintf("%02d\/%02d\/%04d", $month+1, $day, $year+1900);
#my $title = "Distribution of References Spanned by Pattern Size (total refs: $temp, $temp2)";

    #my $imwidth = 2000;
    #my $imheight = 2000;
    #my $graph = GD::Graph::linespoints->new($imwidth, $imheight);

    #  $graph->set(
    #      transparent  => 0,
    #      bgclr    => 'white',
    #      x_label           => 'pattern size',
    #      y_label           => 'number of references',
    #      title             => $title,
    #      x_min_value       => 0,
    #      y_max_value       => int($maxval + $maxval*.1),
    #      x_label_skip      => 10
    #  ) or die $graph->error;

#my @legend_keys = ('input reference','references with at least one spanning read');
#$graph->set_legend(@legend_keys);

    # @data = (
    #    [ @arX0 ],
    #    [ @arY1 ],
    #    [ @arY2 ]
    #  );

    #my $gd = $graph->plot(\@data) or die $graph->error;

    #open(IMG, ">${latex}.span${MIN_SUPPORT_REQUIRED}.psize.png") or die $!;
    #binmode IMG;
    #print IMG $gd->png;
    #close IMG;

    $tall = 0;
    my $tclustered = 0;
    my $tmapped    = 0;
    my $tbest      = 0;
    $tspan1 = 0;

    # 3
    print $distrfh "\n\nArraySize, All, Clustered, PercentageC, MappedFlanks,"
        . " PercentageM, BBB, PercentageB, Span1, PercentageS\n";
    $sth = $dbh->prepare(q{
        SELECT (lastindex - firstindex + 1) AS arraysize, count(distinct rid)
        FROM refdb.fasta_ref_reps
        GROUP BY arraysize
        ORDER BY arraysize ASC})
        or die "Couldn't prepare statement: " . $dbh->errstr();
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH3{ $data[0] } = $data[1];
        }
        else {
            $PHASH3{$LARGEST_ASIZE} += $data[1];
        }
        $i++;
    }

    # 4
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct mainreftab.rid)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.span1>0 GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH4{ $data[0] } = $data[1];
        }
        else {
            $PHASH4{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }

    # bbb and percent mapped (PHASHB and PHASH6)
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize,
    COUNT(reftab.rid), SUM(bbb)
    FROM refdb.fasta_ref_reps reftab JOIN
    (SELECT refid AS rid, MAX(bbb) AS bbb FROM map GROUP BY rid)
    USING (rid)
    GROUP BY arraysize ORDER BY arraysize ASC}
    );
    $sth->execute();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH6{ $data[0] } = $data[1];
            $PHASHB{ $data[0] } = $data[2];
        }
        else {
            $PHASH6{$LARGEST_ASIZE} += $data[1];
            $PHASHB{$LARGEST_ASIZE} += $data[2];
        }

        $i++;
    }

    # percent clustered
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct reftab.rid)
    FROM refdb.fasta_ref_reps reftab INNER JOIN
    (SELECT -repeatid AS rid FROM clusterlnk WHERE repeatid < 0)
    USING(rid)
    GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH5{ $data[0] } = $data[1];
        }
        else {
            $PHASH5{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }

    # print
    $maxval = 10;
    @arX0   = ();
    @arX    = ();
    @arY1   = ();
    @arY2   = ();
    foreach my $key ( sort { $a <=> $b } ( keys %PHASH3 ) ) {
        my $val1 = 0;
        if ( exists $PHASH3{$key} ) { $val1 = $PHASH3{$key}; }
        my $val2 = 0;
        if ( exists $PHASH4{$key} ) { $val2 = $PHASH4{$key}; }
        my $val4 = 0;
        if ( exists $PHASH5{$key} ) { $val4 = $PHASH5{$key}; }
        my $val3 = 0;
        if ( exists $PHASH6{$key} ) { $val3 = $PHASH6{$key}; }
        my $valb = 0;
        if ( exists $PHASHB{$key} ) { $valb = $PHASHB{$key}; }

        if ( $key < $LARGEST_ASIZE ) {
            print $distrfh $key . ", "
                . $val1 . ", "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $valb . ", "
                . int( 100 * $valb / $val1 ) . "%, "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%\n";
            push( @arX0, "$key" );
        }
        else {
            print $distrfh $key . "+, "
                . $val1 . ", "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $valb . ", "
                . int( 100 * $valb / $val1 ) . "%, "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%\n";
            push( @arX0, "$key+" );
        }
        push( @arX,  $key );
        push( @arY1, $val1 );
        push( @arY2, $val2 );

        $tall       += $val1;
        $tclustered += $val4;
        $tmapped    += $val3;
        $tbest      += $valb;
        $tspan1     += $val2;

        $maxval = max( $maxval, $val1 );
    }

    print $distrfh "TOTAL, $tall, $tclustered, "
        . int( 100 * $tclustered / $tall )
        . "\%, $tmapped, "
        . int( 100 * $tmapped / $tall )
        . "\%, $tbest, "
        . int( 100 * $tbest / $tall )
        . "\%, $tspan1, "
        . int( 100 * $tspan1 / $tall ) . "\%\n";

    # copies gained/lost
    my $total = 0;
    print $distrfh
        "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) Copies Gained, Frequency\n";
    $sth = $dbh->prepare(
        qq{SELECT copies, COUNT(*)
        FROM vntr_support
        WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED
        GROUP BY copies ORDER BY copies ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();

    while ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0] . "," . $data[1] . "\n";
        $i++;
        $total += $data[1];
    }
    print $distrfh "TOTAL: $total\n";

    # copies by patsize
    $total = 0;
    print $distrfh
        "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) PatternSize, Copies Gained, Frequency\n";
    $sth = $dbh->prepare(
        qq{SELECT length(pattern) AS patsize, copies, count(*)
        FROM vntr_support INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
        WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED
        GROUP BY patsize, copies
        ORDER BY patsize ASC, copies ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0] . "," . $data[1] . "," . $data[2] . "\n";
        $i++;
        $total += $data[2];
    }
    print $distrfh "TOTAL: $total\n";

    # copies by array size
    $total = 0;
    print $distrfh
        "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) ArraySize, Copies Gained, Frequency\n";
    $sth = $dbh->prepare(
        qq{SELECT (lastindex-firstindex+1) AS arraysize, copies, count(*)
        FROM vntr_support INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
        WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED
        GROUP BY arraysize, copies
        ORDER BY arraysize ASC, copies ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();

    while ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0] . "," . $data[1] . "," . $data[2] . "\n";
        $i++;
        $total += $data[2];
    }
    print $distrfh "TOTAL: $total\n";

    # Counts for TR and allele spanning reads

    # Allele support
    my %AllelesSing   = ();
    my %AllelesIndist = ();
    my %AllelesBBB    = ();

    # my %AllelesDist   = ();

    # Allele support, TRs called VNTRs (not alternate alleles only)
    my %VNTRAllelesSing   = ();
    my %VNTRAllelesIndist = ();
    my %VNTRAllelesBBB    = ();

    # my %VNTRAllelesDist   = ();

    # Ref TR support
    my %SpanningSing   = ();
    my %SpanningIndist = ();
    my %SpanningBBB    = ();

    # my %SpanningDist   = ();

    $sth = $dbh->prepare(
        q{SELECT refid, is_singleton, support_vntr > 0, GROUP_CONCAT(copies),
    GROUP_CONCAT(support)
    FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab
        ON reftab.rid = -vntr_support.refid
    INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
    GROUP BY refid}
    );
    $sth->execute();
    while ( my @data = $sth->fetchrow_array() ) {
        my $tr_support_total = 0;

        # my @copies = split ",", $data[3];
        my @support = split ",", $data[4];
        for my $s ( 0 .. $#support ) {
            $tr_support_total              += $support[$s];
            $AllelesSing{ $support[$s] }   += $data[1];
            $AllelesIndist{ $support[$s] } += !$data[1];
            $AllelesBBB{ $support[$s] }++;
            $VNTRAllelesSing{ $support[$s] }   += ( $data[1]  && $data[2] );
            $VNTRAllelesIndist{ $support[$s] } += ( !$data[1] && $data[2] );
            $VNTRAllelesBBB{ $support[$s] }    += $data[2];
        }
        $SpanningSing{$tr_support_total}   += $data[1];
        $SpanningIndist{$tr_support_total} += !$data[1];
        $SpanningBBB{$tr_support_total}++;
    }

    # SPANNING
    print $distrfh "\n\nSpanning reads per locus\n\nRefClass     ";
    my $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %SpanningBBB ) ) {
        $maxkey = $key;
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nSINGLETON";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanningSing{$key} ) { $val1 = $SpanningSing{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nINDIST      ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanningIndist{$key} ) { $val1 = $SpanningIndist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nBBB         ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanningBBB{$key} ) { $val1 = $SpanningBBB{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    # ALLELES
    print $distrfh "\n\nSpanning reads per allele\n\nRefClass     ";
    $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBB ) ) {
        $maxkey = $key;
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nSINGLETON";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesSing{$key} ) { $val1 = $AllelesSing{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nINDIST      ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesIndist{$key} ) { $val1 = $AllelesIndist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nBBB         ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBB{$key} ) { $val1 = $AllelesBBB{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";


    # ALLELES for VNTRs
    print $distrfh "\n\nSpanning reads per allele for VNTRs\n\nRefClass     ";
    $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %VNTRAllelesBBB ) ) {
        $maxkey = $key;
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nSINGLETON";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $VNTRAllelesSing{$key} ) {
            $val1 = $VNTRAllelesSing{$key};
        }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

# $total = 0;
# print $distrfh "\nDISTING   ";
# for ( my $key = 0; $key <= $maxkey; $key++ ) {
#     my $val1 = 0;
#     if ( exists $VNTRAllelesDist{$key} ) { $val1 = $VNTRAllelesDist{$key}; }
#     print $distrfh "\t$val1";
#     $total += $val1;
# }
# print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nINDIST      ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $VNTRAllelesIndist{$key} ) {
            $val1 = $VNTRAllelesIndist{$key};
        }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nBBB         ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $VNTRAllelesBBB{$key} ) { $val1 = $VNTRAllelesBBB{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total\n";

    # alleles spanning reads for VNTRs by VNTR class
    # TODO Below queries could probably be optimized further
    # (or folded into the query above)

    my %AllelesBBBInf   = ();
    my %AllelesBBBObs   = ();
    my %AllelesBBBTotal = ();

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE support_vntr>0 AND mainreftab.homez_diff=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesBBBInf{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE support_vntr>0 AND (mainreftab.hetez_same=1 OR mainreftab.hetez_diff=1 OR mainreftab.hetez_multi=1)}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesBBBObs{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    print $distrfh
        "\n\nSpanning reads per allele for VNTRs by VNTR class\n\nVNTRClass     ";
    $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBBInf ) ) {
        $maxkey = $key;
        $AllelesBBBTotal{$key} = $AllelesBBBInf{$key};
    }
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBBObs ) ) {
        $maxkey = $key;
        $AllelesBBBTotal{$key} += $AllelesBBBObs{$key};
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nInferred";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBBInf{$key} ) { $val1 = $AllelesBBBInf{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nObserved";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBBObs{$key} ) { $val1 = $AllelesBBBObs{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nTotal";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBBTotal{$key} ) {
            $val1 = $AllelesBBBTotal{$key};
        }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total\n";

    # counts for support 1 alleles
    print $distrfh "\n1A: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE copies!=0 AND mainreftab.has_support=0 AND support_vntr=0 AND support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    print $distrfh "\n1B: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE copies!=0 AND has_support=1 AND support_vntr=0 AND support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    print $distrfh "\n2A: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE homez_diff=1 and support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    print $distrfh "\n2B: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE (hetez_same=1 OR hetez_diff=1 OR hetez_multi=1) AND support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    $dbh->disconnect;
    close($distrfh);

    # images

    #$imwidth = 2000;
    #$imheight = 2000;
    #$graph = GD::Graph::linespoints->new($imwidth, $imheight);

#$title = "Distribution of References Spanned by Array Size (total refs: $temp, $temp2)";

    #  $graph->set(
    #      x_label           => 'array size',
    #      y_label           => 'number of references',
    #      title             => $title,
    #      x_min_value       => 0,
    #      y_max_value       => int($maxval + $maxval*.1),
    #      x_label_skip      => 10
    #  ) or die $graph->error;

#@legend_keys = ('input reference','references with at least one spanning read');
#$graph->set_legend(@legend_keys);

    # @data = (
    #    [ @arX0 ],
    #    [ @arY1 ],
    #    [ @arY2 ]
    #  );

    #$gd = $graph->plot(\@data) or die $graph->error;

    #open(IMG, ">${latex}.span${MIN_SUPPORT_REQUIRED}.asize.png") or die $!;
    #binmode IMG;
    #print IMG $gd->png;
    #close IMG;

}