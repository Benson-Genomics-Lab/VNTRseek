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
#  vntrseek 0 19 --dbsuffix run1 --server orca.bu.edu --nprocesses 8
#    --html_dir /var/www/html/vntrview --fasta_dir /bfdisk/watsontest
#    --output_root /smdisk --tmpdir /tmp &
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
use File::Path "remove_tree";
use File::Temp;

# use POSIX qw(strftime);
use Getopt::Long qw(GetOptionsFromArray);
use List::Util qw(min max);
use lib ( "$FindBin::RealBin/lib", "$FindBin::RealBin/local/lib/perl5" );
use vutil
    qw(get_config validate_config set_statistics get_statistics write_sqlite
       get_dbh get_ref_dbh set_datetime print_config create_blank_file run_redund trim);

my $VERSION = '@VNTRVer@';
my $install_dir = "$FindBin::RealBin"; # where the pipeline is installed
my $run_dir = getcwd();                # where the script gets called from

#my $REFS_REDUND = 374579;
my $REFS_REDUND = 0;    # used to be used in latex file, not anymore

my $HELPSTRING
    = "\nUsage: $0 <start step> <end step> [options]\n\n"
    . "\tThe first step is 0, last step is 19.\n"
    . "\tAt least --DBSUFFIX, --INPUT_DIR, and --REFERENCE must be provided,\n"
    . "\t  or a valid --CONFIG file specifying them.\n"
    . "\t  Use --GEN_CONFIG to generate a default file.\n\n"
    . "\tA config can be provided with command line arguments,\n"
    . "\t  in which case the command line values will take precedence.\n"
    . "\n\nOPTIONS:\n\n"
    . "\t--HELP                        prints this help message\n"
    . "\t--DBSUFFIX                    prefix for database name (such as the name of your analysis)\n"
    . "\t--INPUT_DIR                   input data directory (BAM or plain or gzipped fasta/fastq files)\n"
    . "\t--OUTPUT_ROOT                 output directory (must be writable and executable!)\n"
    . "\t--TMPDIR                      temp (scratch) directory (must be writable!)\n"
    . "\t--SERVER                      server name, used for generating html links\n"
    . "\n"
    . "\t--REFERENCE                   base name of reference files (default set in global config file)\n"
    . "\t--REDO_REFDB                  force reinitialization of the reference set database, 0/1 (default 0)\n"
    . "\t--REFERENCE_INDIST_PRODUCE    generate a file of indistinguishable references, 0/1 (default 0)\n"
    . "\n"
    . "\t--PLOIDY                      sample's ploidy (default 2)\n"
    . "\t--READ_LENGTH                 length of fasta reads (default 150)\n"
    . "\t--IS_PAIRED_READS             data is paired reads, 0/1 (default 1)\n"
    . "\t--STRIP_454_KEYTAGS           for 454 platform, strip leading 'TCAG', 0/1 (default 0)\n"
    . "\t--KEEPPCRDUPS                 whether to find and remove PCR duplicates. (default: 1, duplicates are kept)\n"
    . "\n"
    . "\t--MIN_FLANK_REQUIRED          minimum required flank on both sides for a read TR to be considered (default 10)\n"
    . "\t--MAX_FLANK_CONSIDERED        maximum flank length used in flank alignments, set high to use full flank (default 50)\n"
    . "\t--MIN_SUPPORT_REQUIRED        minimum number of mapped reads which agree on copy number to call an allele (default 2)\n"
    . "\n"
    . "\t--NPROCESSES                  number of processors to use on your system (default 2)\n"
    . "\t--CONFIG                      path to file with any of the above options. Command line option only.\n"
    . "\t--GEN_CONFIG                  set with a path to generate a clean config. Command line option only.\n"
    . "\t--CLEAN                       force reinitialization of the run database. Command line option only.\n"
    . "\t--STATS                       print out a simple table of run statistics. Command line option only.\n"
    . "\t\n\nADDITIONAL USAGE:\n\n"
    . "\t.../vntrseek.pl 100           clear error\n"
    . "\t.../vntrseek.pl 100 N         clear error and set NextRunStep to N (0-19)\n"
    . "\t                                This is only when running on a cluster using the\n"
    . "\t                                advanced cluster script that checks for NextRunStep\n\n"
    . "\t.../vntrseek.pl --REFERENCE </path/to/reference/basename> [--REDO_REFDB]\n"
    . "\t  Initialize a reference set. This only needs to be done once for a reference set.\n"
    . "\t  Supply --REDO_REFDB to force a recreation of the database.\n\n";

# set options based on command line
my %opts = ();
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
    "READ_LENGTH=i",          "CONFIG=s",
    "GEN_CONFIG"
);

# Print help string if that's all they ask
if ($opts{'HELP'} or !@ARGV) { print $HELPSTRING; exit;}

# Generate a copy of the default config to serve as a template.
if ($opts{'GEN_CONFIG'}) {
    if (-d $opts{'GEN_CONFIG'}) { $opts{'GEN_CONFIG'} .= '/example.vs.cnf';}
    die "File $opts{'GEN_CONFIG'} already exists. Will not overwrite."
        if -e $opts{'GEN_CONFIG'};

    system("cp $install_dir/defaults.vs.cnf $opts{'GEN_CONFIG'}");
    print "Generated example configuration file at $opts{'GEN_CONFIG'}.";
    exit;
}


# Load configuration files and combine with CLI args.
# Make sure all the configuration options are safe for running.
%opts = get_config(%opts);
validate_config();

# the output of the pipeline, ex: "/bfdisk/vntr_$DBSUFFIX";
my $output_folder = "$opts{OUTPUT_ROOT}/vntr_$opts{DBSUFFIX}";
mkdir $output_folder unless -e $output_folder;

# Write out a terse config file to output directory.
my $config_file = print_config();


# enter install directory
die("Cannot access install directory ($install_dir).\n")
    unless chdir("$install_dir");

# Initialize a reference set db if --reference was given, but not --dbsuffix
if (!exists $opts{'DBSUFFIX'} && exists $opts{'REFERENCE'} ) {
    print "\nPreparing reference set $opts{'REFERENCE'}...\n";
    my $dbh = get_ref_dbh( $opts{'REFERENCE'}, { redo => $opts{'REDO_REFDB'} } );
    die "Error importing reference set\n" unless $dbh;
    exit;
}

# Print a nicely formatted table with all stats if "STATS" flag is given
if (exists $opts{'STATS'}) {
    # add a check here for database before trying to retrieve it?
    my $dbh = get_dbh();
    my $stat_sth = $dbh->prepare(qq{SELECT * FROM stats WHERE id = 1});
    $stat_sth->execute();
    my @stat_cols = $stat_sth->{'NAME'}->@*;
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

        print sprintf("%-${min_width}s\t%20s", $s, ( $run_stats{$s} // "undef" )) . "\n";
    }
    exit;
}

my $timestart;

# this is where TRF output will go converted to leb36 format
my $read_profiles_folder = "$output_folder/data_out/";

# this is where renumbered and non-cyclicly redundant leb36 files will go
my $read_profiles_folder_clean = "$output_folder/data_out_clean/";

# this is where edges are calculated
my $edges_folder = "$opts{TMPDIR}/vntr_$opts{DBSUFFIX}/edges/";

my $DBNAME     = "VNTRPIPE_$opts{DBSUFFIX}";
my $HTTPSERVER = "$opts{SERVER}/vntrview";


# clustering parameters (only cutoffs, other nonessantial paramters are in run_proclu.pl
my $PROCLU_EXECUTABLE = "psearch.exe";
my $CLUST_PARAMS      = " 88 ";
my $MINPROFSCORE      = .85;

# TRF options
my $MATCH    = 2;
my $MISMATCH = 5;
my $INDEL    = 7;
my $MIN_PERIOD_REQUIRED = 7;    # anything with pattern less then this is discarded

my $TRF_EXECUTABLE = '@TRFBin@';
my $TRF_PARAM
    = "'./$TRF_EXECUTABLE' - $MATCH $MISMATCH $INDEL 80 10 50 2000 -d -h -ngs";

my $TRF2PROCLU_EXE = 'trf2proclu-ngs.exe';
my $TRF2PROCLU_PARAM
    = "'./$TRF2PROCLU_EXE' -m $MATCH -s $MISMATCH -i $INDEL -p $MIN_PERIOD_REQUIRED -l $opts{MIN_FLANK_REQUIRED}";

# verify executables
my @executables = (
    $install_dir, $TRF_EXECUTABLE, $TRF2PROCLU_EXE, $PROCLU_EXECUTABLE, "redund.exe",
    "flankalign.exe", "refflankalign.exe", "pcr_dup.exe", "join_clusters.exe");

for my $exec (@executables) {
    die("'$exec' not executable!") unless (-e $exec and -x _);
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
        print "ClearError: making next step: $to.\n\n";
        my $set_step = min(max( $to, 0), 19 );
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
            print min( $idx + 1, 20 ) . "\n";
            return 0;
        }
    }
    return 0;
}

####################################

# set to 0 to do one step at a time (recommended for test run), 1 to run though all steps
#(THIS IS OTIONAL AS SPECIFYING END STEP IS POSSIBLE)
my $DOALLSTEPS = 0;
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

# pipeline error checking

if ( $STEP != 0 ) {

    # clear error?
    if ( $STEP == 100 ) {
        if ($STEPEND >= 0 && $STEPEND <= 19) { ClearError($STEPEND);}
        else                                 { ClearError();}

        warn "Pipeline error cleared!\n";
        exit;
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

    print "Executing step #$STEP (creating database).\n";

    # This function checks if the output directory exists,
    # and creates it if it doesn't exist.
    write_sqlite();

    if ( !-e "$read_profiles_folder" && !mkdir("$read_profiles_folder") ) {
        warn "Warning: Failed to create output directory!\n";
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

    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 1; }

}

if ( $STEP == 1 ) {
    unlink glob "${read_profiles_folder}/*.clu";

    print "Executing step #$STEP (searching for tandem repeats in reads,"
         . " producing profiles and sorting)...\n";
    my $extra_param  = ( $opts{STRIP_454_KEYTAGS} ) ? '-s' : '';
    my $extra_param2 = ( $opts{IS_PAIRED_READS} )   ? '-r' : '';

    $timestart = time();
    my $exstring = "$install_dir/run_trf_ng.pl"
                   . qq( -t "$TRF_PARAM")
                   . qq( -u "$TRF2PROCLU_PARAM")
                   . " $extra_param $extra_param2"
                   . " -p $opts{NPROCESSES}"
                   . qq( "$opts{INPUT_DIR}" "$read_profiles_folder");
    my $trf_out = qx($exstring);
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 2; }
}

if ( $STEP == 2 ) {

    print "Executing step #$STEP (reassigning IDs to repeats).\n";
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

    print "Done!\n\n";

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
    # END TODO
    print "Executing step #$STEP (eliminating cyclic redundancies).\n";

    if ( !mkdir("$read_profiles_folder_clean") ) {
        warn "\nWarning: Failed to create output directory!\n";
    }

    # Run redund.exe on read TRs
    system("./redund.exe",
        $read_profiles_folder,
        "$read_profiles_folder_clean/allreads.leb36",
        "-i");
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

    print "Setting additional statistics...\n";
    system("./setdbstats.pl",
        "$read_profiles_folder",
        "$read_profiles_folder_clean",
        "$config_file");
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 4; }
}

if ( $STEP == 4 ) {

    unlink glob "${read_profiles_folder_clean}/*.clu";
    unlink glob "${read_profiles_folder_clean}/*.cnf";
    unlink glob "${read_profiles_folder_clean}/*.proclu_log";

    $timestart = time();
    print "Executing step #$STEP (bipartite clustering of TR profiles).\n";

    system("./checkleb36.pl",
        $read_profiles_folder_clean);
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
        print $tmp_filt_file_fh join( " ", @fields ) . "\n";
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

    system("./run_proclu.pl",
        1,
        $read_profiles_folder_clean,
        $reference_folder,
        $CLUST_PARAMS,
        $opts{'NPROCESSES'},
        $PROCLU_EXECUTABLE,
        0,
        $opts{'MAX_FLANK_CONSIDERED'});
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
    print "Removing clusters with no references (clara/pam split).\n";
    opendir( DIR, $read_profiles_folder_clean );
    my @files = grep( /clu$/, readdir(DIR) );
    closedir(DIR);

    foreach my $file (@files) {
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 5; }
}

if ( $STEP == 5 ) {

    $timestart = time();
    print "Executing step #$STEP (joining clusters from different proclu runs on reference ids).\n";

    system("./join_clusters.exe",
        $read_profiles_folder_clean,
        "$read_profiles_folder_clean/all.clusters");
    if ( $? == -1 ) {
        SetError( $STEP, "command failed: $!", -1 );
        die "command failed: $!\n";
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) {
            SetError($STEP,
                "joining clusters from different proclu runs on reference ids failed",
                $rc
            );
            die "command exited with value $rc";
        }
    }

    set_statistics( { TIME_JOINCLUST => time() - $timestart } );
    set_datetime("DATE_JOINCLUST");
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 6; }
}

if ( $STEP == 6 ) {
    print "Step #$STEP is empty.\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 7; }
}

if ( $STEP == 7 ) {
    $timestart = time();
    print "\nStep #$STEP is empty.\n\n";

    set_statistics( { TIME_DB_INSERT_REFS => time() - $timestart } );
    set_datetime("DATE_DB_INSERT_REFS");

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 8; }
}

if ( $STEP == 8 ) {

    print "Executing step #$STEP (inserting READS flanks into database).\n";

    $timestart = time();

    # Get rotindex saved in db
    my $dbh = get_ref_dbh( $opts{REFERENCE}, { redo => $opts{REDO_REFDB} } );
    my ($rotindex_str) = $dbh->selectrow_array('SELECT rotindex FROM files');
    $dbh->disconnect();
    die "Error getting rotindex file. Try rerunning with --redo_refdb option.\n"
        unless $rotindex_str;

    my $reference_folder = File::Temp->newdir();
    open my $tmp_rotindex, ">", "$reference_folder/reference.leb36.rotindex";

    # Need to negate all indices
    $rotindex_str =~ s/(\d+)/-$1/g;
    print $tmp_rotindex $rotindex_str;
    close $tmp_rotindex;

    system("./insert_reads.pl",
        "$read_profiles_folder_clean/all.clusters",
        $read_profiles_folder,
        $opts{'INPUT_DIR'},
        $read_profiles_folder_clean,
        "$reference_folder/reference.leb36.rotindex",
        $opts{'STRIP_454_KEYTAGS'},
        $config_file,
        $opts{'TMPDIR'},
        $opts{'IS_PAIRED_READS'});
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 9; }
}

if ( $STEP == 9 ) {

    $timestart = time();
    print "Executing step #$STEP (outputting flanks inside each cluster).\n";

    system("./run_flankcomp.pl",
        "$read_profiles_folder_clean/allwithdups.clusters",
        $opts{'DBSUFFIX'},
        $config_file,
        $opts{'TMPDIR'},
        "$read_profiles_folder_clean/allwithdups.flanks");
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 10; }
}

if ( $STEP == 10 ) {

    $timestart = time();

    system("rm -Rf ${read_profiles_folder_clean}/out");

    if ( !mkdir("${read_profiles_folder_clean}/out") ) {
        warn "\nWarning: Failed to create output directory!\n";
    }

    print "Executing step #$STEP (aligning ref-read flanks).\n";

# 0 for maxerror, means flankalign will pick maxerror based on individual flanklength
    system("./flankalign.exe",
        "$read_profiles_folder_clean/out",
        "$read_profiles_folder_clean/result",
        "$read_profiles_folder_clean/allwithdups.flanks",
        0,
        $opts{'MAX_FLANK_CONSIDERED'},
        $opts{'NPROCESSES'},
        15);
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 11; }
}

if ( $STEP == 11 ) {
    $timestart = time();
    print "Step #$STEP is empty.\n\n";

    set_statistics( { TIME_MAP_REFFLANKS => time() - $timestart } );
    set_datetime("DATE_MAP_REFFLANKS");

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 12; }
}

if ( $STEP == 12 ) {

    $timestart = time();

    system("rm -Rf ${read_profiles_folder_clean}/result");

    if ( !mkdir("${read_profiles_folder_clean}/result") ) {
        warn "\nWarning: Failed to create output result directory!\n";
    }

    print "Executing step #$STEP (inserting map and rankflank information into database).\n";

    system("./run_rankflankmap.pl",
        "$read_profiles_folder_clean/allwithdups.clusters",
        "$read_profiles_folder_clean/out",
        $opts{'TMPDIR'},
        $config_file);
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 13; }
}

if ( $STEP == 13 ) {

    $timestart = time();
    print "Executing step #$STEP (calculating edges).\n";

    system("./run_edges.pl",
        $opts{'REFERENCE'},
        $edges_folder,
        $config_file,
        $MINPROFSCORE,
        $opts{'NPROCESSES'},
        $PROCLU_EXECUTABLE,
        $opts{'TMPDIR'});
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

    rmdir $edges_folder;

    set_statistics( { TIME_EDGES => time() - $timestart } );
    set_datetime("DATE_EDGES");
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 14; }
}

if ( $STEP == 14 ) {

    $timestart = time();
    print "Executing step #$STEP (generating .index files for pcr_dup).\n";

    system("./extra_index.pl",
        "$read_profiles_folder_clean/best",
        $config_file);
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 15; }
}

if ( $STEP == 15 ) {

    $timestart = time();
    print "Executing step #$STEP (calculating PCR duplicates).\n";
    if ( $opts{KEEPPCRDUPS} ) {
        print "Notice: Not removing detected PCR duplicates\n";
    }

    system("./pcr_dup.pl",
        "$read_profiles_folder_clean/best",
        $read_profiles_folder_clean,
        $opts{'DBSUFFIX'},
        $config_file,
        $opts{'NPROCESSES'},
        $opts{'TMPDIR'},
        $opts{'KEEPPCRDUPS'});
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 16; }
}

if ( $STEP == 16 ) {

    $timestart = time();
    print "Executing step #$STEP (removing mapped duplicates).\n";

    system("./map_dup.pl",
        $config_file,
        $opts{'TMPDIR'},
        "$read_profiles_folder_clean/result/$opts{DBSUFFIX}.map_dup.txt");
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 17; }
}

if ( $STEP == 17 ) {

    $timestart = time();
    print "Executing step #$STEP (computing variability).\n";

    system("./run_variability.pl",
        "$read_profiles_folder_clean/allwithdups.clusters",
        "$read_profiles_folder_clean/out",
        $config_file,
        $opts{'MIN_FLANK_REQUIRED'},
        $opts{'TMPDIR'});
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
    print "Done!\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 18; }
}

if ( $STEP == 18 ) {
    $timestart = time();
    print "Step #$STEP is empty.\n\n";

    set_statistics( { TIME_ASSEMBLYREQ => time() - $timestart } );
    set_datetime("DATE_ASSEMBLYREQ");

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ($DOALLSTEPS)         { $STEP = 19; }
}

if ( $STEP == 19 ) {

    $timestart = time();
    print "Executing step #$STEP (final database update).\n";

    unless ( get_statistics(qw(NUMBER_TRS_IN_READS)) ) {

        # lets do this setdbstats again (sometimes when copying
        # databases steps are omited so this might not have been
        # executed)
        print "Setting additional statistics.\n";
        system("./setdbstats.pl",
            "$read_profiles_folder",
            "$read_profiles_folder_clean",
            "$config_file");
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

    # further output files (vcf)
    system("./updaterefs.pl",
        $read_profiles_folder,
        $read_profiles_folder_clean,
        $opts{'DBSUFFIX'},
        $config_file,
        "${read_profiles_folder_clean}/result/${DBNAME}",
        $VERSION);
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
    ) if ( exists $opts{HTML_DIR}
           && !( -e "$opts{HTML_DIR}/$opts{DBSUFFIX}/result" ) );

    # cleanup
    print "Cleanup Time!\n";

    remove_tree("${read_profiles_folder_clean}/best", {safe => 1});
    remove_tree("${read_profiles_folder_clean}/out", {safe => 1});

    set_statistics( { TIME_REPORTS => time() - $timestart } );
    set_datetime("DATE_REPORTS");

    # Create reduced database
    my $dbfile = "$output_folder/$opts{DBSUFFIX}.db";
    my $dbfile2 = "$output_folder/$opts{DBSUFFIX}_rl$opts{READ_LENGTH}.db";
    print "Finalizing:\n$dbfile\nbeing reduced into\n$dbfile2\n";
    system("sqlite3 $dbfile < reduced_db.sql");
    system("mv temp_reduced.db $dbfile2");

    print "Done!\n\n";
}

print "Finished!";
