#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use POSIX qw(strftime);
use Carp qw(croak carp);
use FindBin;
use File::Basename;
use Try::Tiny;
use Data::Dumper;
use lib "$FindBin::RealBin/lib";

use vutil
    qw(get_config get_dbh get_trunc_query set_statistics get_statistics gen_exec_array_cb vs_db_insert);

my $RECORDS_PER_INFILE_INSERT = 100000;

#use GD::Graph::linespoints;

=head1

updaterefs.pl - Calculates VNTRs from database and writes
final reports and VCF files

=cut

# Argument parsing & set up
my $argc = @ARGV;
die "Usage: updaterefs.pl expects 7 arguments.\n"
    unless $argc >= 7;

my $curdir        = getcwd();
my $readpf        = $ARGV[0]; # read_profiles_folder
my $rpfc          = $ARGV[1]; # read_profiles_folder_clean
my $DBSUFFIX      = $ARGV[2]; # dbname
my $cnf           = $ARGV[3]; # config_file
my $filerep       = $ARGV[4]; # file_representative
my $result_prefix = $ARGV[5]; # latexfile
my $VERSION       = $ARGV[6]; # VERSION

my %run_conf = get_config("CONFIG", $cnf);
my ( $HTTPSERVER, $MIN_SUPPORT_REQUIRED, $TEMPDIR )
    = @run_conf{qw(SERVER MIN_SUPPORT_REQUIRED TMPDIR)};

############################ Procedures ###############################################################

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
    croak "make_vcf_rec requires one argument"
        unless @_ == 1;
    my $supported_tr = shift;

    my $qual = ".";

    my $filter = ( $supported_tr->{is_singleton} == 1 ) ? "PASS" : "SC";

    my $info = sprintf(
        "RC=%.2lf;RPL=%d;RAL=%d;RCP=%s;ALGNURL=http://%s/index.php?db=VNTRPIPE_%s&ref=-%d&isref=1&istab=1&ispng=1&rank=3",
        $supported_tr->{copynum}, length( $supported_tr->{consenuspat} ),
        $supported_tr->{arlen},   $supported_tr->{consenuspat},
        $HTTPSERVER,              $DBSUFFIX,
        $supported_tr->{rid}
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

#open (VCFFILE2, ">", "${latex}.span2.vcf") or die "Can't open for reading ${latex}.vcf.";

####################################
sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

####################################
sub commify {
    my $input = shift;
    carp "Undef input" unless defined $input;
    $input = reverse $input;
    $input =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return reverse $input;
}

####################################
# Takes a boolean as an argument. If the boolean is anything perl considers false,
# then this function will only produce a VCF file for supported VNTRs. If true, all
# supported TRs are included in the file.

sub print_vcf {

    # my $allwithsupport = shift;

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
    my $numsup_sth
        = $dbh->prepare(
        "select count(distinct vntr_support.refid) FROM vntr_support WHERE support >= $MIN_SUPPORT_REQUIRED"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;

    my $numsup = 0;
    $numsup_sth->execute() or die "Cannot execute: " . $numsup_sth->errstr();
    $numsup_sth->bind_columns( \$numsup );
    $numsup_sth->fetch;
    if ( !defined $numsup ) {
        die "Error getting number of supported TRs: " . $dbh->errstr;
    }
    $numsup_sth->finish;

    # Get number of VNTRs supported
    my $numvntrs_sth
        = $dbh->prepare(
        "select count(*) FROM fasta_ref_reps WHERE support_vntr>0")
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
##database=VNTRPIPE_$DBSUFFIX
##databaseurl=http://${HTTPSERVER}/result.php?db=VNTRPIPE_${DBSUFFIX}
##INFO=<ID=RC,Number=1,Type=Float,Description="Reference Copies">
##INFO=<ID=RPL,Number=1,Type=Integer,Description="Reference Pattern Length">
##INFO=<ID=RAL,Number=1,Type=Integer,Description="Reference Tandem Array Length">
##INFO=<ID=RCP,Number=1,Type=String,Description="Reference Consensus Pattern">
##INFO=<ID=ALGNURL,Number=1,Type=String,Description="Alignment URL">
##FILTER=<ID=SC,Description="Reference is Singleton">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SP,Number=R,Type=Integer,Description="Number of Spanning Reads">
##FORMAT=<ID=CGL,Number=R,Type=Integer,Description="Copies Gained or Lost with respect to reference">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$DBSUFFIX
];

    print $spanN_vcffile $vcf_header;
    print $allwithsupport_vcffile $vcf_header;

    # Note: This query simply fetches all representative read TRs
    # for all supported alleles of all supported ref TRs, including
    # the reference allele (if supported).
    # This means ORDER is important. If changes to the query need
    # to be made, pay attention that the reference supporting read
    # is be the FIRST one for each TR. Then the reads for the
    # remaining alleles, in increasing copy number.
    # Because the ordering of items in GROUP_CONCAT in SQLite is
    # not predictable, I use a subquery here instead to ensure
    # the ordering I want.
    my $get_supported_reftrs_query
        = qq{SELECT rid, is_singleton, firstindex, arlen, copynum,
    pattern AS consenuspat, head, REPLACE(UPPER(sequence), " ", "") AS sequence,
    refdir, GROUP_CONCAT(copies) AS cgl, (MAX(sameasref) == 1) AS refdetected,
    GROUP_CONCAT(support) AS support, GROUP_CONCAT(readarray) AS alt_seqs,
    GROUP_CONCAT(readdir) AS readdir, COUNT(*) AS num_alleles
    FROM (SELECT reftab.rid, is_singleton, firstindex,
        (lastindex - firstindex) + 1 AS arlen, reftab.copynum, reftab.pattern,
        reftab.head, sequence, c1.direction AS refdir, copies, sameasref, support,
        REPLACE(UPPER(SUBSTR(dna, first, (last-first)+1)), " ", "") AS readarray,
        c2.direction AS readdir
        FROM main.fasta_ref_reps mainreftab
        INNER JOIN vntr_support ON mainreftab.rid=-vntr_support.refid
        INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid=-vntr_support.refid
        INNER JOIN clusterlnk c1 ON vntr_support.refid=c1.repeatid
        INNER JOIN replnk ON vntr_support.representative=replnk.rid
        INNER JOIN clusterlnk c2 ON c2.repeatid=replnk.rid
        INNER JOIN fasta_reads ON replnk.sid=fasta_reads.sid
        WHERE support >= ${MIN_SUPPORT_REQUIRED}
        ORDER BY reftab.head ASC, reftab.firstindex ASC, sameasref DESC, copies ASC
    ) GROUP BY rid};
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
            die
                "Error: sequence not found in database for ref ($row{rid}) alternate allele(s): ",
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

        $row{sequence}
            = ( $row{sequence} ) ? $row{sequence} : ".";
        $row{support}   = join( ",", @read_support );
        $row{cgl}       = join( ",", @cgl );
        $row{alt}       = join( ",", @rev_seqs );
        $row{gt_string} = join( "/", @gt );

        my $vcf_rec = make_vcf_rec( \%row );

        # Only print record to spanN file if VNTR
        ( $row{num_alleles} > 1 || !$row{refdetected} )
            && $vntr_count++
            && print $spanN_vcffile $vcf_rec;
        print $allwithsupport_vcffile $vcf_rec;
    }

    warn "Warning: Mismatch in VNTR count. Supported vntr count is $numvntrs "
        . "but we counted $vntr_count when producing VCF files. VCF header "
        . "will have a bad VNTR count.\n"
        unless ( $vntr_count == $numvntrs );
    close $spanN_vcffile;
    close $allwithsupport_vcffile;
    $dbh->disconnect;
}

####################################
sub print_distr {
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
        or die
        "\nCan't open for reading ${result_prefix}.span${MIN_SUPPORT_REQUIRED}.txt\n";
    print $distrfh
        "\n\nPatternSize, All, Span1, PercentageS1, Span${MIN_SUPPORT_REQUIRED}, PercentageS${MIN_SUPPORT_REQUIRED}, VntrSpan${MIN_SUPPORT_REQUIRED}, PercentageVS${MIN_SUPPORT_REQUIRED}\n";

    my $dbh = get_dbh( { readonly => 1, userefdb => 1 } );

    # 1 patsize (all)
    my $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM refdb.fasta_ref_reps GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
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
    $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.span1>0
            GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
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
    $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.spanN>0 GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
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
    $sth = $dbh->prepare(
        q{SELECT LENGTH(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE support_vntr>0 GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
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
    warn "TOTAL, $tall, $tspan1, \n";
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
    print $distrfh
        "\n\nArraySize, All, Clustered, PercentageC, MappedFlanks, PercentageM, BBB, PercentageB, Span1, PercentageS\n";
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct rid)
        FROM refdb.fasta_ref_reps
        GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
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

    #     q{SELECT refid, SUM(support), is_singleton
    # FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab
    #     ON reftab.rid = -vntr_support.refid
    # INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
    # GROUP BY refid}
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

# $sth = $dbh->prepare(
#     q{SELECT refid,sum(support)
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
#         INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
#     WHERE is_dist=1 GROUP BY refid}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# my %SpanninDist;

    # while ( my @data = $sth->fetchrow_array() ) {
    #     $SpanninDist{ $data[1] }++;
    #     $i++;
    # }
    # $sth->finish;

# $sth = $dbh->prepare(
#     q{SELECT refid,sum(support)
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
#         INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
#     WHERE is_indist=1 GROUP BY refid}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $SpanningIndist{ $data[1] }++;
#     $i++;
# }
# $sth->finish;

# $sth = $dbh->prepare(
#     q{SELECT refid,sum(support)
#     FROM vntr_support INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
#     GROUP BY refid}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $SpanningBBB{ $data[1] }++;
#     $i++;
# }
# $sth->finish;

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

    # $total = 0;
    # print $distrfh "\nDISTING   ";
    # for ( my $key = 0; $key <= $maxkey; $key++ ) {
    #     my $val1 = 0;
    #     if ( exists $SpanninDist{$key} ) { $val1 = $SpanninDist{$key}; }
    #     print $distrfh "\t$val1";
    #     $total += $val1;
    # }
    # print $distrfh "\t$total";

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

    # # alleles spanning reads

# $sth = $dbh->prepare(
#     q{SELECT refid,copies,support
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
#         INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
#     WHERE is_singleton=1}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $AllelesSing{ $data[2] }++;
#     $i++;
# }
# $sth->finish;

# # $sth = $dbh->prepare(
# #     q{SELECT refid,copies,support
# #     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
# #         INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
# #     WHERE is_dist=1}
# # ) or die "Couldn't prepare statement: " . $dbh->errstr;
# # $sth->execute() or die "Cannot execute: " . $sth->errstr();
# # while ( my @data = $sth->fetchrow_array() ) {
# #     $AllelesDist{ $data[2] }++;
# #     $i++;
# # }
# # $sth->finish;

# $sth = $dbh->prepare(
#     q{SELECT refid,copies,support
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
#         INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
#     WHERE is_indist=1}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $AllelesIndist{ $data[2] }++;
#     $i++;
# }
# $sth->finish;

# $sth = $dbh->prepare(
#     q{SELECT refid,copies,support
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $AllelesBBB{ $data[2] }++;
#     $i++;
# }
# $sth->finish;

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

    # $total = 0;
    # print $distrfh "\nDISTING   ";
    # for ( my $key = 0; $key <= $maxkey; $key++ ) {
    #     my $val1 = 0;
    #     if ( exists $AllelesDist{$key} ) { $val1 = $AllelesDist{$key}; }
    #     print $distrfh "\t$val1";
    #     $total += $val1;
    # }
    # print $distrfh "\t$total";

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

    # alleles spanning reads for VNTRs

    # %AllelesSing   = ();
    # # %AllelesDist   = ();
    # %AllelesIndist = ();
    # %AllelesBBB    = ();

# $sth = $dbh->prepare(
#     q{SELECT refid,copies,support
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
#         INNER JOIN refdb.fasta_ref_reps reftab USING(rid)
#     WHERE  support_vntr>0 AND is_singleton=1}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $AllelesSing{ $data[2] }++;
#     $i++;
# }
# $sth->finish;

# # $sth = $dbh->prepare(
# #     q{SELECT refid,copies,support
# #     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
# #         INNER JOIN refdb.fasta_ref_reps reftab USING(rid)
# #     WHERE support_vntr>0 AND is_dist=1}
# # ) or die "Couldn't prepare statement: " . $dbh->errstr;
# # $sth->execute() or die "Cannot execute: " . $sth->errstr();
# # while ( my @data = $sth->fetchrow_array() ) {
# #     $AllelesDist{ $data[2] }++;
# #     $i++;
# # }
# # $sth->finish;

# $sth = $dbh->prepare(
#     q{SELECT refid,copies,support
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
#         INNER JOIN refdb.fasta_ref_reps reftab USING(rid)
#     WHERE support_vntr>0 AND is_indist=1}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $AllelesIndist{ $data[2] }++;
#     $i++;
# }
# $sth->finish;

# $sth = $dbh->prepare(
#     q{SELECT refid,copies,support
#     FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
#     WHERE support_vntr>0}
# ) or die "Couldn't prepare statement: " . $dbh->errstr;
# $sth->execute() or die "Cannot execute: " . $sth->errstr();
# while ( my @data = $sth->fetchrow_array() ) {
#     $AllelesBBB{ $data[2] }++;
#     $i++;
# }
# $sth->finish;

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

###################
sub print_latex {
    my $ReadTRsSupport = shift;

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
# TIE-OK -- applies to a read which is MAP, but may appear in more than one flank TTT cluster (due to ties or no flanking sequence)
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
        = $stat_hash->{CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS}
        ;    # INITIAL READ-TRs GE7 PC ADDBACK

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
     #my $refTRsMappedSingleton = GetStatistics("NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED");  # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SINGLETON

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
        == $readTRsMappedToIndistinguishable + $readTRsMappedToSingleton )
    {
        print Dumper(
            [   $readTRsMapped,
                $readTRsMappedToIndistinguishable,
                $readTRsMappedToSingleton
            ]
        ) . "\n";
        die
            "Error: mismatch of sum of mapped read TRs by distinguishablility and total read TRs mapped. "
            . "(Expected $readTRsMapped, got "
            . (
            $readTRsMappedToIndistinguishable + $readTRsMappedToSingleton )
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

    unless ( $refTRsMapped
        == $refTRsMappedIndistinguishable + $refTRsMappedSingleton )
    {
        warn Dumper(
            \(  $refTRsMappedSingleton,
                $InvarTRsMappedSingleton,
                $VNTRasSingleton,
                $refTRsMappedIndistinguishable,
                $InvarTRsMappedIndistinguishable,
                $VNTRasIndistinguishable
            )
        );
        die
            "Error: mismatch of sum of mapped ref TRs by distinguishablility and total ref TRs mapped. "
            . "(Expected $refTRsMapped, got '$refTRsMappedIndistinguishable' + '$refTRsMappedSingleton' = '"
            . ( $refTRsMappedIndistinguishable + $refTRsMappedSingleton )
            . "')\n";
    }

    my $refWithSupport
        = $sum_has_support;   # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT
    my $refAsHomozygousSame = $sum_homez_same
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HOMZYGOUS SAME
    my $refAsHomozygousDiff = $sum_homez_diff
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HOMZYGOUS DIFF
    my $refAsHeterozygousSame = $sum_hetez_same
        ;   # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HETEROZYGOUS SAME
    my $refAsHeterozygousDiff = $sum_hetez_diff
        ;   # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HETEROZYGOUS DIFF
    my $refAsMulti = $sum_hetez_multi
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT MULTI

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
            . "\n\\end{center}", $DBSUFFIX
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

    $sth
        = $dbh->prepare(
        "SELECT count(distinct head) FROM map INNER JOIN replnk ON map.readid=replnk.rid INNER JOIN fasta_reads ON fasta_reads.sid=replnk.sid  WHERE bbb=1;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
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

# read representative file
my %DHASH   = ();
my %REPHASH = ();

#open (MYFILE, "$filerep") or die "Can't open for reading $filerep.";
#while (<MYFILE>) {
#
#  if ($_ =~ /(\d+)_(\d+)\W(.+)\W(.+)/i) {
#    $id=$1;
#    $fd=$3;
#    $pt=$4;
#print $id. " " . $fd. " ". $pt . "\n";
#    $DHASH{$id} = $fd;
#    $REPHASH{$id} = $pt;
#  }
#}

set_statistics( { N_MIN_SUPPORT => $MIN_SUPPORT_REQUIRED } );

my $dbh = get_dbh( { userefdb => 1 } )
    or die "Could not connect to database: $DBI::errstr";

#goto AAA;

$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");

# Store temporary tables in memory
$dbh->do("PRAGMA temp_store = 2");

# $dbh->begin_work;
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "main.fasta_ref_reps" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

# $dbh->commit;

# Create a temporary table of all refids from vntr_support
# table, with the data aggregated to make parsing easier.
# Must negate refid to match in JOIN later.
$dbh->do(
    q{CREATE TEMPORARY TABLE temp.vntr_support AS
    SELECT rid,pattern,sameasref,support
    FROM refdb.fasta_ref_reps INNER JOIN (
        SELECT -refid AS rid, GROUP_CONCAT(sameasref) AS sameasref,
        GROUP_CONCAT(support) AS support
        FROM vntr_support GROUP BY refid ORDER BY -refid ASC
    ) USING (rid)}
);

# Avoid DB lock issues with shared database
$dbh->do(q{DETACH DATABASE refdb});

my ($supported_vntr_count)
    = $dbh->selectrow_array(
    q{SELECT COUNT(DISTINCT rid) FROM temp.vntr_support});

# now we can SELECT from reference TR table, using the temp table
# to JOIN.
my $get_supported_reftrs_sth = $dbh->prepare(
    q{SELECT *
    FROM temp.vntr_support}
) or die "Couldn't prepare statement: " . $dbh->errstr;

warn "\nUpdating fasta_ref_reps table...\n";

my $query = qq{INSERT INTO main.fasta_ref_reps
        (rid, alleles_sup, allele_sup_same_as_ref, entropy,
            has_support, span1, spanN, homez_same, homez_diff,
            hetez_same, hetez_diff, hetez_multi, support_vntr,
            support_vntr_span1)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    };

my $update_ref_table_sth = $dbh->prepare($query)
    or die "Couldn't prepare statement: " . $dbh->errstr;

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

    ( $ENV{DEBUG} ) && warn "Saving ref entry: " . Dumper($ref) . "\n";
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
            "Error when inserting entries into our supported reference TRs table. Row dump: "
                . Dumper( \@supported_refTRs ) );
        if ($rows) {
            @supported_refTRs = ();
        }
        else {
            die
                "Something went wrong inserting, but somehow wasn't caught!\n";
        }
    }
}

my $updfromtable = $i;

# Insert last rows:
if (@supported_refTRs) {
    my $cb   = gen_exec_array_cb( \@supported_refTRs );
    my $rows = vs_db_insert( $dbh, $update_ref_table_sth, $cb,
        "Error when inserting entries into our supported reference TRs table. Row dump: "
            . Dumper( \@supported_refTRs ) );
    if ($rows) {
        @supported_refTRs = ();
    }
    else {
        die "Something went wrong inserting, but somehow wasn't caught!\n";
    }
}
if ( $updfromtable != $supported_vntr_count ) {
    die
        "Updated number of entries($updfromtable) not equal to the number of references, aborting!";
}

# updating stats table
print STDERR "Updating stats table...\n";
my $update_stats_sth = $dbh->prepare(
    q{UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED=(
    SELECT COUNT(*)
        FROM (
            SELECT COUNT(*) AS thecount
            FROM clusterlnk
            WHERE repeatid<0 GROUP BY clusterid HAVING thecount=1
        ) f
    )}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$dbh->begin_work;
$update_stats_sth->execute()
    or die "Cannot execute: " . $update_stats_sth->errstr();
$dbh->commit;

$dbh->begin_work;
$dbh->do('CREATE TEMPORARY TABLE t1 (c1 INT PRIMARY KEY NOT NULL);')
    or die "Couldn't do statement: " . $dbh->errstr;
$dbh->commit;

my $sth = $dbh->prepare(
    q{INSERT INTO t1 SELECT urefid
    FROM (
        SELECT COUNT(*) AS thecount,MAX(repeatid) AS urefid
        FROM clusterlnk
        WHERE repeatid<0
        GROUP BY clusterid
        HAVING thecount=1
    ) f}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$dbh->begin_work;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
$dbh->commit;

$sth = $dbh->prepare(
    q{UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER = (
        SELECT COUNT(DISTINCT repeatid)
        FROM clusterlnk 
        WHERE repeatid IN (SELECT c1 FROM t1)
    )}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$dbh->begin_work;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
$dbh->commit;

$sth = $dbh->prepare(
    q{UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED = (
        SELECT COUNT(distinct map.refid)
        FROM clusterlnk INNER JOIN map ON map.refid=-clusterlnk.repeatid
        WHERE repeatid IN (select c1 from t1)
    )}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$dbh->begin_work;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
$dbh->commit;

$sth = $dbh->prepare(
    q{UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED =
        NUMBER_REFS_SINGLE_REF_CLUSTER -
        NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$dbh->begin_work;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
$dbh->commit;

# Update last few stats:
my ( $mapped, $rank, $rankflank, $num_spanN );

($mapped) = $dbh->selectrow_array(q{SELECT COUNT(*) FROM map})
    or die "Couldn't select map count: " . $dbh->errstr;

($rank) = $dbh->selectrow_array(q{SELECT COUNT(*) FROM rank})
    or die "Couldn't select rank count: " . $dbh->errstr;

($rankflank) = $dbh->selectrow_array(q{SELECT COUNT(*) FROM rankflank})
    or die "Couldn't select rankflank count: " . $dbh->errstr;

# update spanN number on stats
($num_spanN) = $dbh->selectrow_array(
    q{SELECT COUNT(*) FROM fasta_ref_reps 
    WHERE support_vntr > 0}
) or die "Couldn't select span N count: " . $dbh->errstr;

$dbh->begin_work;
$dbh->do(
    qq{UPDATE stats SET
        NUMBER_MAPPED=$mapped,
        NUMBER_RANK=$rank,
        NUMBER_RANKFLANK = $rankflank,
        NUMBER_REFS_VNTR_SPAN_N = $num_spanN}
) or die "Couldn't do statement: " . $dbh->errstr;
$dbh->commit;

# set old db settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");
$dbh->do("PRAGMA temp_store = 0");

# Index for getting copies gained/lost
# $dbh->do(q{CREATE INDEX IF NOT EXISTS
#     "idx_vntr_support_copies_support"
#     ON "vntr_support" (`copies`,`support`)}
# );
# Optimize queries
$dbh->do("PRAGMA main.optimize");
$dbh->disconnect;

warn "Producing output LaTeX file...\n";
print_latex($ReadTRsSupport);

warn "Producing distribution file...\n";
print_distr();

# Print VCF
warn "Producing VCF...\n";
print_vcf();

# All with support
# warn "Producing VCF (all with support)...\n";
# print_vcf(1);

warn "Done!\n";
