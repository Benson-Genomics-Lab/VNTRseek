#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(all any);
use DBI;

use FindBin;
use lib "$FindBin::RealBin";
use Product;

use Time::HiRes qw(time tv_interval);

my $t0 = time();

my $n = 35;
my @exclude = (39, 78, 92, 96); # skip ' N \ `

my %xlet_to_code = ("", ""); # the empty case
my %code_to_xlet = ("", "");

# generate all triplets, doublets, singlets
for my $reps (3, 2, 1) {
    my $prod = Product->new($reps, qw(A C G T));
    while (my @set = $prod->next()) {
        my $xlet = join("", @set);

        while (any {$_ == $n} @exclude) { $n++; }
        my $code = chr($n);
        $xlet_to_code{$xlet} = $code;
        $code_to_xlet{$code} = $xlet;
        $n += 1;
    }
}

# alignment encoding token parsing regex
my $decom_re = qr/(\d+)((?:D+)|S|I)([ACGTN]*)/;

# dna reverse complement
my %WCPair = qw(A T C G G C T A);
sub reverse_complement {
    my $out = "";
    for my $c ( reverse split(//, uc $_[0]) ) {
        $out .= defined $WCPair{$c} ? $WCPair{$c}: "N";
    }
    return $out;
}

sub clean_read {
    my $out = uc $_[0];
    $out =~ s/[^ACGT]/N/;
    return $out;
}


# define triplet encode and decode routines
sub triplet_encode {
    # Scan through read and look for non-ACGT characters, which will be replaced with N.
    #   Everything else is replaced by triplet, doublet, or singlet code;
    #   doublet and singlet are used at the end if the sequence length is not a multiple of three
    #   of if a non-ACGT character is encountered.
    #
    # This encoding begins with ! to differentiate from alignmnet encoding.

    my $read = uc $_[0];

    my $encoded_read = "!";
    my $xlet = "";
    my @nucs = qw(A C G T);

    for my $char (split //, $read) {
        if (!any {$_ eq $char} @nucs) {
            $encoded_read .= $xlet_to_code{$xlet} . "N";
            $xlet = "";
            next;
        }
        $xlet .= $char;
        if (length($xlet) == 3) {
            $encoded_read .= $xlet_to_code{$xlet};
            $xlet = "";
        }
    }

    $encoded_read .= $xlet_to_code{$xlet};
    return $encoded_read;
}

sub triplet_decode {
    # Partial reverse of the encoding, since casing and non-ACGT characters are lost.
    # The result is all uppercase, and any non-ACGT in the original are now N.
    my $encoded_read = $_[0];

    # decode
    my $decoded_read = "";
    for my $char (split //, (substr $encoded_read, 1)) {
        $decoded_read .= $char eq "N" ? "N" : $code_to_xlet{$char};
    }

    return $decoded_read;
}


# define align encode and decode routines
sub align_encode {
    # Aligns the read to the reference and encodes it based on the edit distance.
    #   Takes a clip of the reference based on the read flanks.
    # Returns the encoded read and the start position in the reference.

    my ($read, $flankleft, $sequence, $flankright,
        $first, $last, $direction) = @_;
    my $readlength = length $read;

    $read = uc $read;
    my $RC = "";
    if ($direction eq '"') { # if the read should be reverse complement
        $read = reverse_complement($read);
        $RC = "R";
        my $t = $readlength - $last;
        $last = $readlength - $first;
        $first = $t;
    }

    my $clipl = $first >= 10 ? 20+$first : 2*$first;

    my $clipr = ($readlength - $last + 1);
    $clipr = $clipr >= 10 ? 20+$clipr : 2*$clipr;

    my $left_clip = substr $flankleft, -$clipl;
    my $right_clip = substr $flankright, 0, $clipr;

    my $ref_clip = $left_clip . $sequence . $right_clip;
    $ref_clip = uc $ref_clip;

    my $output = qx(./edlib-align $read $ref_clip);
    if (!$output) {
        print $!;
        return "", 0;
    }
    my ($start, $ref_align, $read_align) = split /\n/, $output;
    $start = $start + length($flankleft) - $clipl; # adjust for clip

    # encode read from align strings
    my @readcs = split //, $read_align;
    my @refcs  = split //, $ref_align;
    my ($readc, $refc);

    my $encoded_read = "";
    my $prev = "";
    my $i = 0;
    for $i (0 .. length($read_align) - 1) {
        $readc = $readcs[$i];
        $refc = $refcs[$i];

        if ($readc eq $refc) {
            $prev = ""; next;
        }
        if ($readc eq "-") {
            if ($prev ne "D") { $encoded_read .= "$i"; }
            $encoded_read .= "D";
            $prev = "D";
        }
        elsif ($refc ne "-") {
            if ($prev ne "S") { $encoded_read .= "${i}S"; }
            $encoded_read .= $readc;
            $prev = "S";
        }
        else {
            if ($prev ne "I") { $encoded_read .= "${i}I"; }
            $encoded_read .= $readc;
            $prev = "I";
        }
    }

    if ($encoded_read eq "") { $encoded_read = "="; }
    $encoded_read = $RC . $encoded_read; # reverse complement flag

    return $encoded_read, $start;
}

sub align_decode {

    my ($encoded_read, $flankleft, $sequence, $flankright,
        $start, $readlength, $direction) = @_;

    my $ref = $flankleft . $sequence . $flankright;
    my $decoded_read = substr $ref, $start;

    # decode
    my $RC = substr($encoded_read, 0, 1) eq "R" ? 1 : 0;
    if ($RC) { $encoded_read = substr $encoded_read, 1; }

    if (substr($encoded_read, 0, 1) eq "=") {
        $decoded_read = substr $decoded_read, 0, $readlength;
        $decoded_read = uc $decoded_read;
        return $RC ? reverse_complement($decoded_read) : $decoded_read;
    }

    my $position = -1;
    my $kind = '';
    my $chars = '';
    while ($encoded_read =~ /$decom_re/g) {
        my $position = int($1);
        my $kind = $2;
        my $chars = $3;

        if ($kind eq 'S') {
            $decoded_read = substr($decoded_read, 0, $position)
                . $chars . substr($decoded_read, $position + length($chars));
        }
        elsif ($kind eq 'I') {
            $decoded_read = substr($decoded_read, 0, $position)
                . $chars . substr($decoded_read, $position);
        }
        else {
            $decoded_read = substr($decoded_read, 0, $position)
                . join('', ('-') x length($kind))
                . substr($decoded_read, $position + length($kind));
        }
    }

    $decoded_read =~ s/-//g;
    $decoded_read = substr($decoded_read, 0, $readlength);
    $decoded_read = uc $decoded_read;

    return $RC ? reverse_complement($decoded_read) : $decoded_read;
}


# get name of db
my $dbf   = $ARGV[0]; # database with fasta_reads.reads to compress
my $ref   = $ARGV[1]; # reference database
my $tmp   = $ARGV[2]; # temp folder, see RUNNING.md#computing-clusters
my $procs = $ARGV[3]; # number of child processes to use

my $BATCH = 10000;    # batch size for children, this ↑ => RAM ↑

# connect to db
my $dbh = DBI->connect("DBI:SQLite:dbname=$dbf");
$dbh->do("ATTACH DATABASE ? AS refdb", {}, $ref);

# do big query
# Explanation:
#  The scores table has a lot of duplicates of readids,
#  causing many rows in the intermediate table.
#
#  In addition, the left join from fasta_reads ---to ensure all
#  reads are captured--- disrupts the index used for the group by.
#
#  Together, this creates a prohibitively slow query.
#  To avoid this, we reduce the duplication in scores first
#  with the uncorrelated table t that's grouped by replnk.rid = readid.
#  Since one replnk.rid only has one sid, this greatly simplifies
#  the left join and doesn't disturb the indexing.
#
# QUERY PLAN
# |--MATERIALIZE T
# |  |--SCAN scores
# |  |--SEARCH replnk USING INTEGER PRIMARY KEY (rowid=?)
# |  |--SEARCH clusterlnk USING INTEGER PRIMARY KEY (rowid=?)
# |  |--USE TEMP B-TREE FOR GROUP BY
# |  `--USE TEMP B-TREE FOR count(DISTINCT)
# |--SCAN fasta_reads
# |--SEARCH T USING AUTOMATIC COVERING INDEX (sid=?)
# `--SEARCH refdb.fasta_ref_reps USING INTEGER PRIMARY KEY (rowid=?)
my $data_sth = $dbh->prepare(q{
    SELECT TOTAL(refn), sid, fasta_reads.head, dna,
        flankleft, sequence, flankright, first, last, direction
    FROM fasta_reads
      LEFT JOIN (
          SELECT sid, COUNT(DISTINCT refid) as refn, refid, first, last, direction
          FROM replnk
            JOIN scores     ON replnk.rid = scores.readid
            JOIN clusterlnk ON replnk.rid = clusterlnk.repeatid
          GROUP BY replnk.rid
      ) AS t USING(sid)
      LEFT JOIN refdb.fasta_ref_reps ON t.refid = refdb.fasta_ref_reps.rid
    GROUP BY sid
    ORDER BY sid
});

$data_sth->execute();


my $coconut = -1;
my %p = ();
my $fh;

my $encoded_read;
my $encoded_start;

sub fork_compress {

    my @queue = ();
    my @data;
    for (my $k = 0; $k < $BATCH; $k++) {
        @data = $data_sth->fetchrow_array();
        if (@data == 0) { last;}
        push @queue, [@data];
    }
    if (@queue == 0) { return 0;}

    $coconut++;
    defined(my $pid = fork) or die "Can't fork!";
    if ($pid != 0) { return $pid;}  # ← parent, child ↓

    my ($mult, $sid, $head, $read, $flankleft, $flankright,
        $sequence, $first, $last, $direction);

    my $decoded_read;
    my $original;

    open $fh, ">", "$scratch/${coconut}_compressed_reads.txt";
    for my $k (@queue) {
        ($mult, $sid, $head, $read, $flankleft, $sequence,
             $flankright, $first, $last, $direction) = @$k;

        $original = clean_read($read);
        my $read_length = length($read);

        # 1000
        if (($mult != 1) or (length($sequence) > 1000 - 2*$read_length)) {
            $encoded_read = triplet_encode($read);
            $encoded_start = -1;

            $decoded_read = triplet_decode($encoded_read);
        }
        else {
            ($encoded_read, $encoded_start) = align_encode($read,
                $flankleft, $sequence, $flankright,
                $first, $last, $direction);

            $decoded_read = align_decode($encoded_read,
                $flankleft, $sequence, $flankright,
                $encoded_start, $read_length, $direction);
        }

        ($original eq $decoded_read) or die "Read $sid original not equal to decoded read.\n";
        print $fh join("\t", $sid, $head, $encoded_start, $read_length, $encoded_read) . "\n";
    }
    exit 0; # child does not return
}

for (my $k = 0; $k < $procs; $k++) { $p{fork_compress()} = $coconut;}

while ((my $pid = wait) != -1) {

    # check return value
    my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
    if ($core) {
        warn "process $pid dumped core\n";
        exit(1000);
    }
    elsif ($sig == 9) {
        warn "process $pid was murdered!\n";
        exit(1001);
    }
    elsif ($rc != 0) {
        warn "process $pid has returned $rc!\n";
        exit($rc);
    }

    # one instance has finished processing, start a new one
    delete $p{$pid};
    $p{fork_compress()} = $coconut;
}
$data_sth->finish();
$dbh->disconnect();

# Renew connection after fork
$dbh = DBI->connect("DBI:SQLite:dbname=$dbf");
$dbh->do("PRAGMA synchronous = OFF");

$dbh->do("DROP TABLE IF EXISTS encoded_reads");
$dbh->do(q{
    CREATE TABLE encoded_reads (
        sid integer NOT NULL,
        head integer NOT NULL,
        start integer DEFAULT NULL,
        read_length integer DEFAULT NULL,
        compressed_read varchar(250) DEFAULT NULL,
        PRIMARY KEY(sid) )
});


# Insert compressed results
my $ins_sth = $dbh->prepare("INSERT INTO encoded_reads VALUES (?, ?, ?, ?, ?)");
my @insert;

$dbh->begin_work();
for my $k (0..$coconut) {
    open $fh, "<", "$wf/${k}_compressed_reads.txt";
    while (<$fh>) {
        chomp;
        @insert = split /\t/;
        $ins_sth->execute(@insert) or die "Insert Failed for " . join(", ", @insert);
    }
    close $fh;
    unlink "$wf/${k}_compressed_reads.txt";
}
$dbh->commit();


# Final checks
my ($countO) = $dbh->selectrow_array("select count(*) from fasta_reads");
my ($countN) = $dbh->selectrow_array("select count(*) from encoded_reads");

print "$countN stored / $countO original\n";

if ($countO != $countN) {
    warn "Not all reads compressed; keeping original fasta_reads table.\n";
}
else {
    print "Reads successfully converted; removing original fasta_reads table.\n";
    # remove original table when done
    $dbh->do("DROP TABLE IF EXISTS fasta_reads");
}

$dbh->do("PRAGMA synchronous = ON");
$dbh->disconnect();
