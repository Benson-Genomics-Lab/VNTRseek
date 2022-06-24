#!/usr/bin/env perl

# MASTER SCRIPT TO RUN THE TR VARIANT SEARCH PIPELINE
#
# DO NOT USE SPACES IN PATHS AND DO NOT USE DOTS (.) OR HYPHENS (-) IN DBSUFFIX
#
# command line usage example:
#  vntrseek N K --dbsuffix dbsuffix
#       where N is the start step to execute (0 is the first step)
#       and K is the end step (19 is the last step)
#
# example:
#  vntrseek 0 19 --dbsuffix run1 --server orca.bu.edu --nprocesses 8 --html_dir /var/www/html/vntrview --fasta_dir /bfdisk/watsontest --output_root /smdisk --tmpdir /tmp &
#
# special commands:
#  vntrseek 100 --dbsuffix dbsuffix
#       clear error
#  vntrseek 99 --dbsuffix dbsuffix
#       return next step that needs to be run (this can be
#       used for multi/single processor execution flow control used with
#       advanced cluster script)
#  vntrseek 100 N --dbsuffix dbsuffix
#       clear error and set NextRunStep to N (for advanced cluster script)
#
# IMPORTANT: for correct execution, please add these lines to
# [mysqld] section of my.cnf file and restart mysql process:
#
# innodb_buffer_pool_size=1G
# innodb_additional_mem_pool_size=20M

use strict;
use warnings;
use Cwd;
use v5.24;
use FindBin;
use File::Basename;
use File::Copy;
use File::Path qw(remove_tree);
use File::Temp;

# use POSIX qw(strftime);
use Getopt::Long qw(GetOptionsFromArray);
use List::Util qw(min max);
use lib ( "$FindBin::RealBin/lib", "$FindBin::RealBin/local/lib/perl5" );
use vutil
    qw(get_config set_config set_statistics get_statistics write_sqlite get_dbh get_ref_dbh set_datetime print_config create_blank_file run_redund trim);

# VNTRSEEK Version
my $VERSION = '@VNTRVer@';

# this is where the pipeline is installed
my $install_dir = "$FindBin::RealBin";

# this is the working directory of the run when calling the script
my $run_dir = getcwd();

my $DOALLSTEPS = 0
    ; # set to 0 to do one step at a time (recommended for test run), 1 to run though all steps (THIS IS OTIONAL AS SPECIFYING END STEP IS POSSIBLE)

# TRF options
my $TRF_EXECUTABLE = '@TRFBin@';
my $TRF_EXE        = "$install_dir/$TRF_EXECUTABLE";

my $MATCH    = 2;
my $MISMATCH = 5;
my $INDEL    = 7;
my $MIN_PERIOD_REQUIRED
    = 7;    # anything with pattern less then this is discarded

#my $REFS_REDUND = 374579;
my $REFS_REDUND = 0;    # used to be used in latex file, not anymore

# enter install directory
if ( !chdir("$install_dir") ) {
    { die("Install directory does not exist!\n"); }
}

my $HELPSTRING
    = "\nUsage: $0 --DBSUFFIX <run_name> <start step> <end step> To tell the master script what step to execute. The first step is 0, last step is 19."
    . "\n\nOPTIONS:\n\n"
    . "\t--HELP                        prints this help message\n"
    . "\t--DBSUFFIX                    suffix for database name (this is the name of your analysis)\n"
    . "\t--NPROCESSES                  number of processors on your system\n"
    . "\t--READ_LENGTH                 length of fasta reads\n"
    . "\t--MIN_FLANK_REQUIRED          minimum required flank on both sides for a read TR to be considered (default 10)\n"
    . "\t--MAX_FLANK_CONSIDERED        maximum flank length used in flank alignments, set to big number to use full flank (default 50)\n"
    . "\t--MIN_SUPPORT_REQUIRED        minimum number of mapped reads which agree on copy number to call an allele (default 2)\n"
    . "\t--KEEPPCRDUPS                 if set to 0, PCR duplicates will be looked for and removed. (default: PCR duplicates are not removed)\n"
    . "\t--SERVER                      server name, used for html generating links\n"
    . "\t--STRIP_454_KEYTAGS           for 454 platform, strip leading 'TCAG', 0/1 (default 0)\n"
    . "\t--IS_PAIRED_READS             data is paired reads, 0/1 (default 1)\n"
    . "\t--PLOIDY                      sample's ploidy (default 2)\n"
    . "\t--REDO_REFDB                  force reinitialization of the reference set database, 0/1 (default 0)\n"
    . "\t--CLEAN                       force reinitialization of the run database. Command line option only.\n"
    . "\t--HTML_DIR                    html directory (optional, must be writable and executable!)\n"
    . "\t--INPUT_DIR                   input data directory (BAM or plain or gzipped fasta/fastq files)\n"
    . "\t--OUTPUT_ROOT                 output directory (must be writable and executable!)\n"
    . "\t--TMPDIR                      temp (scratch) directory (must be writable!)\n"
    . "\t--REFERENCE                   base name of reference files (default set in global config file)\n"
    . "\t--REFERENCE_INDIST_PRODUCE    generate a file of indistinguishable references, 0/1 (default 0)\n"
    . "\t\n\nADDITIONAL USAGE:\n\n"
    . "\t$0 100                       clear error\n"
    . "\t$0 100 N                     clear error and set NextRunStep to N (0-19, this is only when running on a cluster using the advanced cluster script that checks for NextRunStep)\n"
    . "\t$0 --REFERENCE </path/to/reference/basename> [--REDO_REFDB] to initialize a reference set.\n"
    . "This only needs to be done once for a reference set. Supply --REDO_REFDB to force a recreation of the database.\n\n";

my %opts = ();

# set options based on command line
GetOptions(
    \%opts,                   "HELP",
    "NPROCESSES=i",           "MIN_FLANK_REQUIRED=i",
    "MAX_FLANK_CONSIDERED=i", "MIN_SUPPORT_REQUIRED=i",
    "KEEPPCRDUPS=i",          "DBSUFFIX=s",
    "SERVER=s",               "STRIP_454_KEYTAGS=i",
    "IS_PAIRED_READS=i",      "PLOIDY=i",
    "REDO_REFDB",             "INPUT_DIR=s",
    "OUTPUT_ROOT=s",          "TMPDIR=s",
    "REFERENCE=s",            "REFERENCE_INDIST_PRODUCE=i",
    "CLEAN",                  "STATS",
    "READ_LENGTH=i"
);

# print help if asked
if ( $opts{"HELP"} ) {
    print $HELPSTRING;
    exit 0;
}

if ( !exists $opts{"DBSUFFIX"} ) {
    # Initialize a reference set db if --reference was given, but not --dbsuffix
    if ( exists $opts{REFERENCE} ) {
        print "\nPreparing reference set $opts{REFERENCE}...\n";
        my $dbh
            = get_ref_dbh( $opts{REFERENCE}, { redo => $opts{REDO_REFDB} } );
        die "Error importing reference set\n" unless $dbh;
        exit;
    }

    die("Please set run name (--DBSUFFIX) variable using command line.\n");
}

# locating local config
# Convert any input for DBSUFFIX to lower case.
$opts{"DBSUFFIX"} = lc( $opts{"DBSUFFIX"} );
my $config_file = "$run_dir/" . $opts{'DBSUFFIX'} . ".vs.cnf";

# set variables from vs.cnf
my %conf_vars = get_config( $opts{'DBSUFFIX'}, $run_dir );
for my $k ( keys %conf_vars ) {
    $opts{$k} = $conf_vars{$k}
        if ( !exists $opts{$k} );
}

# Print a nicely formatted table with all stats if "STATS" flag is given
if ( exists $opts{STATS} ) {
    my $dbh = get_dbh();
    my $stat_sth = $dbh->prepare(qq{SELECT * FROM stats WHERE id = 1});
    $stat_sth->execute;
    my @stat_cols = $stat_sth->{NAME}->@*;
    my $min_width = max(map { length($_) } @stat_cols);
    my %run_stats = $stat_sth->fetchrow_hashref()->%*;

    for my $s (@stat_cols) {
        if ($s =~ /^TIME_/ && defined $run_stats{$s}) {
            my ($d, $h, $m);
            $d = int($run_stats{$s} / (24*3600));
            $run_stats{$s} %= (24*3600);
            $h = int($run_stats{$s} / 3600);
            $run_stats{$s} %= 3600;
            $m = int($run_stats{$s} / 60);
            $run_stats{$s} %= 60;
            $run_stats{$s} = sprintf "%02d:%02d:%02d", $h, $m, $run_stats{$s};
            $run_stats{$s} = "$d days " . $run_stats{$s} if ($d);
        }

        say sprintf("%-${min_width}s\t%20s", $s, ( $run_stats{$s} // "undef" ));
    }
    exit;
}

# Die if no starting step is supplied
die "$HELPSTRING\n" unless @ARGV;

# Send config back to vutil lib and validate, with CLI options taking precedent.
set_config(%opts);

# write out the config file
#   no..... why would you do this except to initialize a local copy?
#   you don't want this if one already exists.
# print_config($run_dir);

# print Dumper(\%opts);
# print Dumper(@ARGV);

my $STEP = shift;
my $STEPEND = -1;
if (   (@ARGV)
    && ( $ARGV[0] =~ /^\d+?$/ )
    && $ARGV[0] >= 0
    && $ARGV[0] <= 99 )
{
    $STEPEND = shift;
    if ( $STEPEND > $STEP ) { $DOALLSTEPS = 1; }
}

my $timestart;
my $DBNAME     = "VNTRPIPE_$opts{DBSUFFIX}";
my $HTTPSERVER = "$opts{SERVER}/vntrview";

# clustering parameters (only cutoffs, other nonessantial paramters are in run_proclu.pl
my $PROCLU_EXECUTABLE = "psearch.exe";
my $CLUST_PARAMS      = " 88 ";
my $MINPROFSCORE      = .85;

my $output_folder
    = "$opts{OUTPUT_ROOT}/vntr_$opts{DBSUFFIX}";    # DO NOT CHANGE, this

# will be created at the output root by the pipeline, ex: "/bfdisk/vntr_$DBSUFFIX";

# this is where TRF output will go converted to leb36 format
my $read_profiles_folder = "$output_folder/data_out/";

# this is where renumbered and cyclicly removed reundancy leb36 files will go
my $read_profiles_folder_clean = "$output_folder/data_out_clean/";

# this is where edges are calculated
my $edges_folder = "$opts{TMPDIR}/vntr_$opts{DBSUFFIX}/edges/";

my $TRF_PARAM
    = "'$TRF_EXE' - $MATCH $MISMATCH $INDEL 80 10 50 2000 -d -h -ngs";
my $TRF2PROCLU_EXE = 'trf2proclu-ngs.exe';
my $TRF2PROCLU_PARAM
    = "'./$TRF2PROCLU_EXE' -m $MATCH -s $MISMATCH -i $INDEL -p $MIN_PERIOD_REQUIRED -l $opts{MIN_FLANK_REQUIRED}";

die("Please set doallsteps (DOALLSTEPS) variable.\n") unless ($DOALLSTEPS >= 0)

# verify executables
my @executables = [
    $install_dir, $TRF_EXECUTABLE, $TRF2PROCLU_EXE, "redund.exe", 
    "flankalign.exe", "refflankalign.exe", "pcr_dup.exe", "join_clusters.exe"];

for my $exec (@executables) {
    die("'$exec' not found!") unless (-e $exec);
    die("'$exec' not executable!") unless (-x $exec);
}

################################################################
sub SetError {
    my $argc = @_;
    die "SetError expects 3 parameters, passed $argc.\n" unless ($argc >= 3);

    my $VALUE1 = $_[0];
    my $VALUE2 = $_[1];
    my $VALUE3 = $_[2];

    set_statistics(
        {   ERROR_STEP => $VALUE1,
            ERROR_DESC => $VALUE2,
            ERROR_CODE => $VALUE3,
        }
    );

    return 0;
}

sub ClearError {
    my @stats = qw(
        DATE_MYSQLCREATE
        DATE_TRF
        DATE_RENUMB
        DATE_REDUND
        DATE_PROCLU
        DATE_JOINCLUST
        DATE_JOINCLUST
        DATE_DB_INSERT_REFS
        DATE_DB_INSERT_READS
        DATE_WRITE_FLANKS
        DATE_MAP_FLANKS
        DATE_MAP_REFFLANKS
        DATE_MAP_INSERT
        DATE_EDGES
        DATE_INDEX_PCR
        DATE_PCR_DUP
        DATE_MAP_DUP
        DATE_VNTR_PREDICT
        DATE_ASSEMBLYREQ
        DATE_REPORTS
    );

    if (@_) {
        my $to = int shift;
        print "\nClearError: making next step: $to.\n\n";
        my $set_step = min( $to, 19 );
        $set_step = max( $to, 0 );
        my %clear_stats = map { $_ => undef } @stats[ $set_step .. 19 ];
        set_statistics( \%clear_stats );

        # Fill in the previous step with some dummy date/time if needed
        if ( get_statistics( $stats[ $set_step - 1 ] ) eq "" ) {
            set_datetime( $stats[ $set_step - 1 ] );
        }
    }

    set_statistics(
        {   ERROR_STEP => 0,
            ERROR_DESC => "",
            ERROR_CODE => 0,
        }
    );
    return 0;
}

sub GetError {
    my $err_hash = get_statistics(qw(ERROR_STEP ERROR_DESC ERROR_CODE));
    return $err_hash;
}

sub GetNextStep {
    # note: this will return 0 regardless if step 0 has been completed or not
    my @stats = qw(
        DATE_MYSQLCREATE
        DATE_TRF
        DATE_RENUMB
        DATE_REDUND
        DATE_PROCLU
        DATE_JOINCLUST
        DATE_JOINCLUST
        DATE_DB_INSERT_REFS
        DATE_DB_INSERT_READS
        DATE_WRITE_FLANKS
        DATE_MAP_FLANKS
        DATE_MAP_REFFLANKS
        DATE_MAP_INSERT
        DATE_EDGES
        DATE_INDEX_PCR
        DATE_PCR_DUP
        DATE_MAP_DUP
        DATE_VNTR_PREDICT
        DATE_ASSEMBLYREQ
        DATE_REPORTS
    );

    my %stat_hash = %{ get_statistics(@stats) };
    my @vals      = @stat_hash{@stats};

    if ( $ENV{DEBUG} ) {
        use Data::Dumper;
        warn "GetNextStep vals: " . Dumper( \@vals ) . "\n";
    }

    for ( my $idx = @vals - 1; $idx >= 0; --$idx ) {
        if ( defined $vals[$idx] ) {
            say min( $idx + 1, 20 );
            return 0;
        }
    }
    return 0;
}

####################################

# pipeline error checking

if ( $STEP != 0 ) {

    # clear error?
    if ( $STEP == 100 ) {

        if   ( $STEPEND >= 0 && $STEPEND <= 19 ) { ClearError($STEPEND); }
        else                                     { ClearError(); }

        warn "Pipeline error cleared!\n";
        exit 0;
    }

    my $err_hash = GetError();
    if ( $err_hash->{ERROR_STEP} != 0 ) {
        die "\nPipeline error detected at step $err_hash->{ERROR_STEP} "
            . "(code: $err_hash->{ERROR_CODE}, description: "
            . "'$err_hash->{ERROR_DESC}'). Call this program with step 100 to clear error.\n";
    }

    # get next step
    if ( $STEP == 99 ) {
        exit GetNextStep();
    }

}

####################################

if ( $STEP == 0 ) {

    $timestart = time();

    print STDERR "Executing step #$STEP (creating database)...\n";

    # This function checks if the output directory exists,
    # and creates it if it doesn't exist.
    write_sqlite();

    if ( !-e "$read_profiles_folder" && !mkdir("$read_profiles_folder") ) {
        warn "\nWarning: Failed to create output directory!\n";
    }

    set_statistics(
        {   MAP_ROOT             => $install_dir,
            N_MIN_SUPPORT        => $opts{MIN_SUPPORT_REQUIRED},
            MIN_FLANK_REQUIRED   => $opts{MIN_FLANK_REQUIRED},
            MAX_FLANK_CONSIDERED => $opts{MAX_FLANK_CONSIDERED},
            TIME_MYSQLCREATE     => time() - $timestart,
        }
    );

    set_datetime("DATE_MYSQLCREATE");

    warn "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 1; }

}

if ( $STEP == 1 ) {
    unlink glob "${read_profiles_folder}/*.clu";

    warn
        "\n\nExecuting step #$STEP (searching for tandem repeats in reads, producing profiles and sorting)...\n";
    my $extra_param  = ( $opts{STRIP_454_KEYTAGS} ) ? '-s' : '';
    my $extra_param2 = ( $opts{IS_PAIRED_READS} )   ? '-r' : '';

    $timestart = time();
    my $trf_out = qx(
        $install_dir/run_trf_ng.pl -t "$TRF_PARAM" -u "$TRF2PROCLU_PARAM" $extra_param $extra_param2 -p $opts{NPROCESSES} "$opts{INPUT_DIR}" "$read_profiles_folder"
    );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }

#elsif ($? & 127) {
#    my $rc = ( $? & 127 );
#    SetError( $STEP, "Died with signal $rc while running TRF and TRF2PROCLU", $rc );
#    die sprintf("died with signal %d, %s coredump\n",
#    $rc,  ($? & 128) ? 'with' : 'without');
#}
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "running TRF and TRF2PROCLU", $rc );
            die "command exited with value $rc";
        }
    }

    my %trf_res = (
        reads             => 0,
        num_trs_ge7       => 0,
        num_trs           => 0,
        num_reads_trs_ge7 => 0,
        num_reads_trs     => 0,
    );

    while ( $trf_out =~ /(\w+):(\d+),?/mg ) {
        next unless exists $trf_res{$1};
        $trf_res{$1} = $2;
    }

    set_statistics(
        {   PARAM_TRF                => $TRF_PARAM,
            FOLDER_FASTA             => $opts{INPUT_DIR},
            FOLDER_PROFILES          => $read_profiles_folder,
            TIME_TRF                 => time() - $timestart,
            NUMBER_READS             => $trf_res{reads},
            NUMBER_TRS_IN_READS_GE7  => $trf_res{num_trs_ge7},
            NUMBER_TRS_IN_READS      => $trf_res{num_trs},
            NUMBER_READS_WITHTRS_GE7 => $trf_res{num_reads_trs_ge7},
            NUMBER_READS_WITHTRS     => $trf_res{num_reads_trs},
        }
    );
    set_datetime("DATE_TRF");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 2; }
}

if ( $STEP == 2 ) {

    print STDERR "\n\nExecuting step #$STEP (reassigning IDs to repeats)...";
    $timestart = time();

    system("./renumber.pl $read_profiles_folder");
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "calling renumber.pl on reads profiles folder",
                $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics(
        {   FILE_REFERENCE_LEB => $opts{REFERENCE} . ".leb36",
            FILE_REFERENCE_SEQ => $opts{REFERENCE} . ".seq",
            TIME_RENUMB        => time() - $timestart,
        }
    );
    set_datetime("DATE_RENUMB");

    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 3; }
}

if ( $STEP == 3 ) {

    if ( -d $read_profiles_folder_clean ) {
        remove_tree( $read_profiles_folder_clean,
            { safe => 1, error => \my $err } );
        if ( $err && @$err ) {
            my $errstr = "Error removing clean read profiles directory: "
                . $err->[0];
            SetError( $STEP, $errstr, -1 );
            die $errstr, "\n";
        }
    }

    $timestart = time();

    # TODO Remove this and/or replace with something that one-time
    # generates an leb36 file with the right ref set seqs.
    say STDERR
        "\n\nExecuting step #$STEP (eliminating cyclic redundancies)...";

    if ( !mkdir("$read_profiles_folder_clean") ) {
        warn "\nWarning: Failed to create output directory!\n";
    }

    # END TODO

    # Run redund.exe on read TRs
    system(
        "./redund.exe $read_profiles_folder $read_profiles_folder_clean/allreads.leb36 -i"
    );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP,
                "calling redund.exe on read profiles folder failed", $rc );
            die "command exited with value $rc";
        }
    }

    # if ( !chdir("$read_profiles_folder_clean") ) {
    #     SetError( $STEP, "could not enter clean profiles directory", -1 );
    #     die("read profiles dir does not exist!");
    # }

    # opendir( my $dirh, $read_profiles_folder_clean );
    # my @files = readdir($dirh);
    # closedir($dirh);
    # foreach my $file (@files) {
    #     if ( $file =~ /^(.*)\.(\d+)$/ ) {
    #         unless (move($file, "$2.$1")) {
    #             my $errstr = "Error moving redund output file $file: $!";
    #             SetError($STEP, $errstr, -1);
    #             die $errstr, "\n";
    #         }
    #     }
    # }

    # if ( !chdir("$install_dir") ) {
    #     SetError( $STEP, "could not enter install directory", -1 );
    #     die("install dir does not exist!\n");
    # }

    warn "setting additional statistics...\n";
    system(
        "./setdbstats.pl $read_profiles_folder $read_profiles_folder_clean $opts{DBSUFFIX} $run_dir"
    );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "calling setdbstats.pl failed", $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics(
        {   FOLDER_PROFILES_CLEAN => $read_profiles_folder_clean,
            TIME_REDUND           => time() - $timestart,
        }
    );
    set_datetime("DATE_REDUND");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 4; }
}

if ( $STEP == 4 ) {

    my $exstring;

    unlink glob "${read_profiles_folder_clean}/*.clu";

    # system($exstring);

    unlink glob "${read_profiles_folder_clean}/*.cnf";

    # system($exstring);

    unlink glob "${read_profiles_folder_clean}/*.proclu_log";

    # system($exstring);

    $timestart = time();
    print STDERR
        "\n\nExecuting step #$STEP (performing bipartite clustering of tandem repeats profiles)...";

    $exstring = "./checkleb36.pl $read_profiles_folder_clean";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "checking leb36 files", $rc );
            die "command exited with value $rc";
        }
    }

# 0 for maxerror, means psearch will pick maxerror based on individual flanklength
# TODO Maybe rewrite run_proclu and psearch to use SQLite db.
# psearch output can still be files but this step will then collect them and
# insert into the db. (Maybe merge in next step's `join_clusters.exe`)
# For now, extract the file we need using the database
    my $dbh = get_ref_dbh( $opts{REFERENCE}, { redo => $opts{REDO_REFDB} } );

    # Filtered file, ordered by min representation as given by redund
    # Must negate rids for refset
    my $get_ref_leb36 = q{SELECT -rid, length(pattern) AS patsize,
    printf("%.2f", copynum), proflen, proflenrc, profile, profilerc, nA, nC, nG, nT,
    printf("%s|%s", upper(substr(flankleft, -60)), upper(substr(flankright, 0, 61))) AS flanks
    FROM fasta_ref_reps JOIN ref_profiles USING (rid)
        JOIN minreporder USING (rid)
    ORDER BY minreporder.idx ASC};
    my $get_filtered_set_profiles_sth = $dbh->prepare($get_ref_leb36);
    my $reference_folder              = File::Temp->newdir();

    $get_filtered_set_profiles_sth->execute;
    open my $tmp_filt_file_fh, ">", "$reference_folder/reference.leb36";
    while ( my @fields = $get_filtered_set_profiles_sth->fetchrow_array ) {
        say $tmp_filt_file_fh join( " ", @fields );
    }
    close $tmp_filt_file_fh;

    # Get rotindex saved in db
    my ($rotindex_str) = $dbh->selectrow_array(
        q{SELECT rotindex
        FROM files}
    );
    die
        "Error getting rotindex file from filtered set. Try rerunning with --redo option\n"
        unless ($rotindex_str);
    $dbh->disconnect();

    # Need to negate all indices
    open my $tmp_rotindex, ">", "$reference_folder/reference.leb36.rotindex";
    $rotindex_str =~ s/(\d+)/-$1/g;
    print $tmp_rotindex $rotindex_str;
    close $tmp_rotindex;

    $exstring
        = "./run_proclu.pl 1 $read_profiles_folder_clean $reference_folder \"$CLUST_PARAMS\" $opts{NPROCESSES} '$PROCLU_EXECUTABLE' "
        . 0
        . " $opts{MAX_FLANK_CONSIDERED}";

    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError(
                $STEP,
                "performing bipartite clustering of tandem repeats profiles failed",
                $rc
            );
            die "command exited with value $rc";
        }
    }

    # remove clusters with no references
    print STDERR
        "\nRemoving clusters with no references (clara/pam split)...\n";
    opendir( DIR, $read_profiles_folder_clean );
    my @files = grep( /clu$/, readdir(DIR) );
    closedir(DIR);

    foreach my $file (@files) {
        print $file. "\n";
        if ( open( my $cluf, "<", "$read_profiles_folder_clean/$file" ) ) {
            open( my $clufout, ">",
                "$read_profiles_folder_clean/$file.clean" );
            while (<$cluf>) {
                if (/-/) { print $clufout $_; }
            }
            close($cluf);
            close($clufout);
            system(
                "mv $read_profiles_folder_clean/$file.clean $read_profiles_folder_clean/$file"
            );
            if ( $? == -1 ) {
                SetError( $STEP, "command failed: $!", -1 );
                die "command failed: $!\n";
            }
            else {
                my $rc = ( $? >> 8 );
                if ( 0 != $rc ) {
                    SetError( $STEP,
                        "renaming files in clean directory failed", $rc );
                    die "command exited with value $rc";
                }
            }
        }
    }

    set_statistics(
        {   PARAM_PROCLU => $CLUST_PARAMS,
            TIME_PROCLU  => time() - $timestart,
        }
    );

    set_datetime("DATE_PROCLU");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 5; }
}

if ( $STEP == 5 ) {

    $timestart = time();
    print STDERR
        "\n\nExecuting step #$STEP (joining clusters from different proclu runs on reference ids)...";
    my $exstring
        = "./join_clusters.exe $read_profiles_folder_clean $read_profiles_folder_clean/all.clusters";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError(
                $STEP,
                "joining clusters from different proclu runs on reference ids failed",
                $rc
            );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_JOINCLUST => time() - $timestart } );
    set_datetime("DATE_JOINCLUST");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 6; }
}

if ( $STEP == 6 ) {

    print STDERR "\n\nSTEP #$STEP IS EMPTY!";
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 7; }
}

if ( $STEP == 7 ) {

    $timestart = time();
    my $exstring;
    warn "Executing step #$STEP (this step is EMPTY!)...";

    set_statistics( { TIME_DB_INSERT_REFS => time() - $timestart } );
    set_datetime("DATE_DB_INSERT_REFS");

    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 8; }
}

if ( $STEP == 8 ) {

    print STDERR
        "\n\nExecuting step #$STEP (inserting READS flanks into database)...";

    $timestart = time();
    my $extra_param = ( $opts{STRIP_454_KEYTAGS} ) ? '1' : '0';

    # Get rotindex saved in db
    my $dbh = get_ref_dbh( $opts{REFERENCE}, { redo => $opts{REDO_REFDB} } );
    my $get_rotindex = q{SELECT rotindex FROM files};
    my ($rotindex_str) = $dbh->selectrow_array($get_rotindex);
    $dbh->disconnect();
    die
        "Error getting rotindex file. Try rerunning with --redo_refdb option.\n"
        unless ($rotindex_str);

    my $reference_folder = File::Temp->newdir();
    open my $tmp_rotindex, ">", "$reference_folder/reference.leb36.rotindex";

    # Need to negate all indices
    $rotindex_str =~ s/(\d+)/-$1/g;
    print $tmp_rotindex $rotindex_str;
    close $tmp_rotindex;

    my $exstring
        = qq{./insert_reads.pl $read_profiles_folder_clean/all.clusters "$read_profiles_folder"  "$opts{INPUT_DIR}" "$read_profiles_folder_clean" "$reference_folder/reference.leb36.rotindex" $extra_param $opts{DBSUFFIX} "$run_dir" $opts{TMPDIR} $opts{IS_PAIRED_READS}};
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "inserting READS flanks into database failed",
                $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_DB_INSERT_READS => time() - $timestart } );
    set_datetime("DATE_DB_INSERT_READS");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 9; }
}

if ( $STEP == 9 ) {

    $timestart = time();
    print STDERR
        "\n\nExecuting step #$STEP (outputting flanks inside each cluster)...";
    my $exstring
        = "./run_flankcomp.pl $read_profiles_folder_clean/allwithdups.clusters $opts{DBSUFFIX} $run_dir $opts{TMPDIR} > $read_profiles_folder_clean/allwithdups.flanks";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "outputting flanks inside each cluster failed",
                $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_WRITE_FLANKS => time() - $timestart } );
    set_datetime("DATE_WRITE_FLANKS");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 10; }
}

if ( $STEP == 10 ) {

    $timestart = time();

    system("rm -Rf ${read_profiles_folder_clean}/out");

    if ( !mkdir("${read_profiles_folder_clean}/out") ) {
        warn "\nWarning: Failed to create output directory!\n";
    }

    print STDERR "\n\nExecuting step #$STEP (aligning ref-read flanks)...";

#my $exstring = "./flankalign.exe $read_profiles_folder_clean/out $read_profiles_folder_clean/result $read_profiles_folder_clean/allwithdups.flanks " . min( 8, int(0.4 * $opts{MIN_FLANK_REQUIRED} + .01)) . " $opts{MAX_FLANK_CONSIDERED} $opts{NPROCESSES} 15";

# 0 for maxerror, means flankalign will pick maxerror based on individual flanklength
    my $exstring
        = "./flankalign.exe $read_profiles_folder_clean/out $read_profiles_folder_clean/result $read_profiles_folder_clean/allwithdups.flanks "
        . 0
        . " $opts{MAX_FLANK_CONSIDERED} $opts{NPROCESSES} 15";

    print STDERR "$exstring\n";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "aligning ref-read flanks failed", $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_MAP_FLANKS => time() - $timestart } );
    set_datetime("DATE_MAP_FLANKS");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 11; }
}

if ( $STEP == 11 ) {

    $timestart = time();
    my $exstring;

    # print STDERR "\n\nExecuting step #$STEP (aligning ref-ref flanks)...";
    print STDERR "\n\nExecuting step #$STEP (This step does nothing!)...";

  # TODO Don't use this. Change so that the new produce_indist script is used.
  # Will require new param
  # if ( $opts{REFERENCE_INDIST_PRODUCE} ) {
  #     print STDERR "\n\n(generating indist file)...";

    #     # enter result dir
    #     if ( !chdir("${read_profiles_folder_clean}/result") ) {
    #         { die("result directory does not exist!"); }
    #     }

#     # copy reference file to result dir
#     # TODO change so that referece file does not need to be in install dir
#     system("cp ${install_dir}/$opts{REFERENCE_FILE} .");
#     if ( $? == -1 ) {
#         SetError( $STEP, "command failed: $!", -1 );
#         die "command failed: $!\n";
#     }
#     else {
#         my $rc = ( $? >> 8 );
#         if ( 0 != $rc ) {
#             SetError(
#                 $STEP,
#                 "copying reference leb36 profile file into reference profiles folder",
#                 $rc
#             );
#             die "command exited with value $rc";
#         }
#     }

    #     $exstring
    #         = "${install_dir}/$PROCLU_EXECUTABLE" . " "
    #         . "${reference_folder}/reference.leb36" . " "
    #         . "$opts{REFERENCE_FILE}" . " "
    #         . "${install_dir}/eucledian.dst" . " "
    #         . $CLUST_PARAMS . " "
    #         . 5
    #         . " 0  -r 50 2> "
    #         . "ref_to_ref.proclu_log";
    #     print "\nrunning: $exstring\n";
    #     system($exstring);
    #     if ( $? == -1 ) {
    #         SetError( $STEP, "command failed: $!", -1 );
    #         die "command failed: $!\n";
    #     }
    #     else {
    #         my $rc = ( $? >> 8 );
    #         if ( 0 != $rc ) {
    #             SetError( $STEP, "aligning ref-ref flanks failed", $rc );
    #             die "command exited with value $rc";
    #         }
    #     }

    #     my $refclusfile = "$opts{REFERENCE_FILE}.clu";
    #     system("mv $refclusfile ${install_dir}/");
    #     if ( $? == -1 ) {
    #         SetError( $STEP, "command failed: $!", -1 );
    #         die "command failed: $!\n";
    #     }
    #     else {
    #         my $rc = ( $? >> 8 );
    #         if ( 0 != $rc ) {
    #             SetError( $STEP, "aligning ref-ref flanks failed", $rc );
    #             die "command exited with value $rc";
    #         }
    #     }
    #     $refclusfile = "${install_dir}/$refclusfile";

    #     # create final indist file
    #     my $indist_withpath = "${install_dir}/$opts{REFERENCE_INDIST}";
    #     open my $tofile, ">", "$indist_withpath" or die $!;
    #     open my $file,   "<", "$refclusfile"     or die $!;
    #     while (<$file>) {
    #         my @values = split( ' ', $_ );
    #         my $repcount = @values - 1;

    #         my $i = 0;
    #         foreach my $val (@values) {

    #             $i++;
    #             $val = trim($val);

    #             if ( my $ch = ( $val =~ m/(\-*\d+)([\-\+])$/g ) ) {
    #                 if ( $1 < 0 && $repcount > 2 ) {
    #                     print $tofile $1 . "\n";
    #                 }
    #             }

    #         }

    #     }
    #     close($file);
    #     close($tofile);

#     # remove working files
#     system(
#         "rm -f ${read_profiles_folder_clean}/result/$opts{REFERENCE_FILE} -f"
#     );
#     system(
#         "rm -f ${read_profiles_folder_clean}/result/ref_to_ref.proclu_log -f"
#     );

    #     # go back to install dir
    #     if ( !chdir("$install_dir") ) {
    #         { die("Install directory does not exist!"); }
    #     }

  #   # update DB
  #   # TODO Replace this with a script which simply sets the flank comparison
  #   # parameters (apparently constant at max_errors = 5 and trim_to = 50, is
  #   # this ever used??)
  #     print STDERR "\n\n(updating database with dist/indist info)...";
  #     $exstring = "./update_indist.pl";
  #     my @exargs = ( qw(-r -k5 -t50 -d), $opts{DBSUFFIX}, "-u", $MSDIR );
  #     system( $exstring, @exargs );
  #     if ( $? == -1 ) {
  #         SetError( $STEP, "command failed: $!", -1 );
  #         die "command failed: $!\n";
  #     }
  #     else {
  #         my $rc = ( $? >> 8 );
  #         if ( 0 != $rc ) {
  #             SetError( $STEP,
  #                 "updating database with dist/undist info failed", $rc );
  #             die "command exited with value $rc";
  #         }
  #     }
  # }

    set_statistics( { TIME_MAP_REFFLANKS => time() - $timestart } );
    set_datetime("DATE_MAP_REFFLANKS");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 12; }
}

if ( $STEP == 12 ) {

    $timestart = time();

    system("rm -Rf ${read_profiles_folder_clean}/result");

    if ( !mkdir("${read_profiles_folder_clean}/result") ) {
        warn "\nWarning: Failed to create output result directory!\n";
    }

    if ( !chdir("$install_dir") ) {
        { die("Install directory does not exist!"); }
    }

    print STDERR
        "\n\nExecuting step #$STEP (inserting map and rankflank information into database.)";
    my $exstring
        = "./run_rankflankmap.pl $read_profiles_folder_clean/allwithdups.clusters $read_profiles_folder_clean/out $opts{TMPDIR} $opts{DBSUFFIX} $run_dir";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError(
                $STEP,
                "inserting map and rankflank information into database failed",
                $rc
            );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_MAP_INSERT => time() - $timestart } );
    set_datetime("DATE_MAP_INSERT");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 13; }
}

if ( $STEP == 13 ) {

    $timestart = time();
    print STDERR "\n\nExecuting step #$STEP (calculating edges)...";
    my $exstring = "./run_edges.pl";
    my @exargs   = (
        "$opts{REFERENCE}",   "$edges_folder",
        "$opts{DBSUFFIX}",    "$run_dir",
        "$MINPROFSCORE",      "$opts{NPROCESSES}",
        "$PROCLU_EXECUTABLE", "$opts{TMPDIR}"
    );
    system( $exstring, @exargs );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "calculating edges failed", $rc );
            die "command exited with value $rc";
        }
    }

    # $exstring = qq(find "${edges_folder}" -type f -delete);
    system(qq(find "${edges_folder}" -type f -delete));

    set_statistics( { TIME_EDGES => time() - $timestart } );
    set_datetime("DATE_EDGES");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 14; }
}

if ( $STEP == 14 ) {

    $timestart = time();
    print STDERR
        "\n\nExecuting step #$STEP (generating .index files for pcr_dup)...";

    my $exstring = "./extra_index.pl";
    my @exargs   = (
        "$read_profiles_folder_clean/best", "$opts{DBSUFFIX}",
        "$run_dir",                         "$opts{TMPDIR}"
    );
    system( $exstring, @exargs );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "generating .index files for pcr_dup failed",
                $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_INDEX_PCR => time() - $timestart } );
    set_datetime("DATE_INDEX_PCR");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 15; }
}

if ( $STEP == 15 ) {

    $timestart = time();
    print STDERR
        "\n\nExecuting step #$STEP (calculating PCR duplicates) ...\n";
    if ( $opts{KEEPPCRDUPS} ) {
        warn "\nWarning: Not removing detected PCR duplicates\n";
    }
    my $exstring = "./pcr_dup.pl";
    my @exargs   = (
        "$read_profiles_folder_clean/best", "$read_profiles_folder_clean",
        "$opts{DBSUFFIX}",                  "$run_dir",
        "$opts{NPROCESSES}",                "$opts{TMPDIR}",
        "$opts{KEEPPCRDUPS}"
    );
    warn "exargs: " . join( " ", @exargs ) . "\n";
    system( $exstring, @exargs );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "calculating PCR duplicates failed", $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_PCR_DUP => time() - $timestart } );
    set_datetime("DATE_PCR_DUP");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 16; }
}

if ( $STEP == 16 ) {

    $timestart = time();
    print STDERR "\n\nExecuting step #$STEP (removing mapped duplicates) ...";

    my $exstring
        = "./map_dup.pl $opts{DBSUFFIX} $run_dir $opts{TMPDIR} > $read_profiles_folder_clean/result/$opts{DBSUFFIX}.map_dup.txt";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "calculating mapped duplicates failed", $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_MAP_DUP => time() - $timestart } );
    set_datetime("DATE_MAP_DUP");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 17; }
}

if ( $STEP == 17 ) {

    $timestart = time();
    print STDERR "\n\nExecuting step #$STEP (computing variability)...";
    my $exstring
        = "./run_variability.pl $read_profiles_folder_clean/allwithdups.clusters $read_profiles_folder_clean/out $opts{DBSUFFIX} $run_dir $opts{MIN_FLANK_REQUIRED} $opts{TMPDIR}";
    system($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "computing variability failed", $rc );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_VNTR_PREDICT => time() - $timestart } );
    set_datetime("DATE_VNTR_PREDICT");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 18; }
}

if ( $STEP == 18 ) {

    $timestart = time();

    print STDERR "\n\nSTEP #$STEP IS EMPTY!";

#print STDERR "\n\nExecuting step #$STEP (computing variability - assembly required)...";
#my $exstring = "./run_assemblyreq.pl $read_profiles_folder_clean/allwithdups.clusters $read_profiles_folder_clean/out $opts{DBSUFFIX} $MSDIR $opts{MIN_FLANK_REQUIRED} $opts{TMPDIR}";
#system($exstring);
#if ( $? == -1 ) { SetError($STEP,"command failed: $!",-1); die "command failed: $!\n"; }
#else {
#  my $rc = ($? >> 8);
#  if ( 0 != $rc ) { SetError($STEP,"computing variability (assembly required) failed",$rc); die "command exited with value $rc"; }
#}

    set_statistics( { TIME_ASSEMBLYREQ => time() - $timestart } );
    set_datetime("DATE_ASSEMBLYREQ");
    print STDERR "done!\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 19; }
}

if ( $STEP == 19 ) {

    $timestart = time();
    print STDERR "\n\nExecuting step #$STEP (final database update)...";

    unless ( get_statistics(qw(NUMBER_TRS_IN_READS)) ) {

        # lets do this setdbstats again (sometimes when copying
        # databases steps are omited so this might not have been
        # executed)
        print STDERR "setting additional statistics...\n";
        system(
            "./setdbstats.pl",             "$read_profiles_folder",
            "$read_profiles_folder_clean", "$opts{DBSUFFIX}",
            "$run_dir"
        );
        if ( $? == -1 ) {
            SetError( $STEP, "command failed: $!", -1 );
            die "command failed: $!\n";
        }
        else {
            my $rc = ( $? >> 8 );
            if ( 0 != $rc ) {
                SetError( $STEP, "calling setdbstats.pl failed", $rc );
                die "command exited with value $rc";
            }
        }
    }

    # distribution, image and latex files
    my $exstring = "perl";
    my @exargs   = (
        "./updaterefs.pl",
        "$read_profiles_folder",
        "$read_profiles_folder_clean",
        "$opts{DBSUFFIX}",
        "$run_dir",
        "$read_profiles_folder_clean/out/representatives.txt",
        "${read_profiles_folder_clean}/result/${DBNAME}",
        "$VERSION"
    );
    system( @exargs );
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError( $STEP, "final database update", $rc );
            die "command exited with value $rc";
        }
    }

    # create symlink to html dir so can be browsed from internet
    symlink(
        "${read_profiles_folder_clean}/result",
        "$opts{HTML_DIR}/$opts{DBSUFFIX}"
        )
        if ( exists $opts{HTML_DIR}
        && !( -e "$opts{HTML_DIR}/$opts{DBSUFFIX}/result" ) );

    # cleanup
    print STDERR "Cleanup...\n";

    $exstring = qq(find "${read_profiles_folder_clean}/best" -type f -delete);
    system($exstring);
    $exstring
        = qq(find "${read_profiles_folder_clean}/edges" -type f -delete);
    system($exstring);
    $exstring = qq(find "${read_profiles_folder_clean}/out" -type f -delete);
    system($exstring);

    set_statistics( { TIME_REPORTS => time() - $timestart } );
    set_datetime("DATE_REPORTS");

    # Create reduced database
    my $opb = "$opts{OUTPUT_ROOT}/vntr_$opts{DBSUFFIX}";
    my $dbfile = "$opb/$opts{DBSUFFIX}.db";
    my $dbfile2 = "$opb/$opts{DBSUFFIX}_rl$opts{READ_LENGTH}.db";
    print STDOUT "Finalizing:\n$dbfile\nbeing reduced into\n$dbfile2\n";
    system("nice --20 sqlite3 $dbfile < reduced_db.sql > $opb/finalization.out 2>&1");
    system("mv temp_reduced.db $dbfile2");
    
    print STDERR "done!\n";
}

print STDERR "\n\nFinished!\n\n";
