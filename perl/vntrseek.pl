#!/usr/bin/env perl

# MASTER SCRIPT TO RUN THE TR VARIANT SEARCH PIPELINE
#
# DO NOT USE SPACES IN PATHS AND DO NOT USE INVALID FILE NAME CHARACTERS IN RUN_NAME (e.g. ? ! /)
#
# command line usage example:
#  vntrseek N K --RUN_NAME RUN_NAME
#       where N is the start step to execute (0 is the first step)
#       and K is the end step (20 is the last step)
#
# example:
#  vntrseek 0 20 --RUN_NAME run1 --SERVER orca.bu.edu --NPROCESSES 8
#     --INPUT_DIR /bfdisk/watsontest --OUTPUT_DIR /smdisk &
#
# special commands:
#  vntrseek 100 --RUN_NAME RUN_NAME
#       clear error
#  vntrseek 99 --RUN_NAME RUN_NAME
#       return next step that needs to be run (this can be
#       used for multi/single processor execution flow control
#       used with advanced cluster script)
#  vntrseek 100 N --RUN_NAME RUN_NAME
#       clear error and set NextRunStep to N (for advanced cluster script)
#

use strict;
use warnings;
use Cwd;
use v5.24;
use File::Basename;
use File::Copy;
use File::Path "remove_tree";
use File::Temp;

use POSIX "strftime";
use Getopt::Long "GetOptionsFromArray";
use List::Util qw(min max);
use DBI;
use DBD::SQLite;
use FindBin;
use lib ( "$FindBin::RealBin/lib", "$FindBin::RealBin/local/lib/perl5" );
use vutil
    qw(get_config validate_config print_config get_dbh get_ref_dbh
    set_statistics get_statistics set_datetime write_sqlite);

my $VERSION = '@VNTRVer@';
my $install_dir = "$FindBin::RealBin"; # where the pipeline is installed
my $run_dir = getcwd();                # where the script gets called from

#my $REFS_REDUND = 374579;
my $REFS_REDUND = 0;    # used to be used in latex file, not anymore

my $HELPSTRING
    = "\nUsage: $0 <start step> <end step> [options]\n\n"
    . "\tThe first step is 0, last step is 20.\n"
    . "\tAt least --RUN_NAME, --INPUT_DIR, and --REFERENCE must be provided,\n"
    . "\t  or a valid --CONFIG file specifying them.\n"
    . "\t  Use --GEN_CONFIG to generate a default file.\n\n"
    . "\tA config can be provided with command line arguments,\n"
    . "\t  in which case the command line values will take precedence.\n"
    . "\n\nOPTIONS:\n\n"
    . "\t--HELP                        prints this help message\n"
    . "\t--RUN_NAME                    prefix for database name (such as the name of your analysis)\n"
    . "\t--INPUT_DIR                   input data directory (BAM or plain or gzipped fasta/fastq files)\n"
    . "\t--OUTPUT_DIR                  output directory (must be writable and executable!)\n"
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
    . "\t--CLEANUP                     clear intermediate files if run successful. Command line option only.\n"
    . "\t--STATS                       print out a simple table of run statistics. Command line option only.\n"
    . "\n\nADDITIONAL USAGE:\n\n"
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
    \%opts,         "CONFIG=s",
    "RUN_NAME=s",   "DBSUFFIX=s",
    "INPUT_DIR=s",  "FASTA_DIR=s",
    "OUTPUT_DIR=s", "OUTPUT_ROOT=s",
    "SERVER=s",
    "REFERENCE=s",   "PLOIDY=i",
    "READ_LENGTH=i", "IS_PAIRED_READS=i",
    "STRIP_454_KEYTAGS=i",
    "MIN_FLANK_REQUIRED=i",   "MAX_FLANK_CONSIDERED=i",
    "MIN_SUPPORT_REQUIRED=i", "KEEPPCRDUPS=i",
    "NPROCESSES=i",
    "REDO_REFDB",   "REFERENCE_INDIST_PRODUCE=i",
    "CLEANUP",
    "HELP", "STATS", "GEN_CONFIG:s"
);

# Print help string if that's all they ask
if ($opts{'HELP'} or !@ARGV) { print $HELPSTRING; exit;}

# Generate a copy of the default config to serve as a template.
if ($opts{'GEN_CONFIG'}) {
    if (-d $opts{'GEN_CONFIG'}) { $opts{'GEN_CONFIG'} .= '/example.vs.cnf';}
    die "File $opts{'GEN_CONFIG'} already exists. Will not overwrite."
        if -e $opts{'GEN_CONFIG'};
 
    if ($opts{'GEN_CONFIG'} == 1) {
        $opts{'GEN_CONFIG'} = 'example.vs.cnf';
    }
    system("cp $install_dir/defaults.vs.cnf $opts{'GEN_CONFIG'}");
    print "Generated example configuration file at $opts{'GEN_CONFIG'}.";
    exit;
}


# Load configuration files and combine with CLI args.
# Make sure all the configuration options are safe for running.
%opts = get_config(%opts);
validate_config();

# the output of the pipeline, ex: "/bfdisk/vntr_$RUN_NAME";
my $output_folder = "$opts{OUTPUT_DIR}/vntr_$opts{RUN_NAME}";
die "Cannot access output folder at $output_folder. Check permissions."
    unless (-r -w -e $output_folder or mkdir $output_folder);

# Write out a terse config file to output directory.
my $config_file = print_config();


# enter install directory
die("Cannot access install directory ($install_dir).\n")
    unless chdir("$install_dir");

# Initialize a reference set db if --reference was given, but not --RUN_NAME
if (!exists $opts{'RUN_NAME'} && exists $opts{'REFERENCE'} ) {
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
my $trff = "$output_folder/data_out/";

# this is where renumbered and non-cyclicly redundant leb36 files will go
my $processedf = "$output_folder/data_out_clean/";

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
    die("'$exec' not executable!") unless (-x -e $exec);
}

my $STEP = shift;
die "Start step must be a positive integer.\n" unless ($STEP =~ /^\d+$/);
my $STEPEND = @ARGV ? shift: -1;
die "End step must be an integer.\n" unless ($STEPEND =~ /^-?\d+$/);

die "Invalid steps given.\n" unless ($STEP == 100 or $STEPEND >= $STEP);
my $CONTINUOUS = ( $STEPEND > $STEP );

################################################################
sub FlagError {
    my $state = shift . " failed";
    my %stats = (ERROR_STEP => $STEP, ERROR_DESC => $state);
    if ( $? == -1 ) {
        $stats{'ERROR_DESC'} .= ": $!";
        $stats{'ERROR_CODE'} = -1;
    }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { $stats{'ERROR_CODE'} = $rc;}
    }
    if (defined $stats{'ERROR_CODE'}) {
        set_statistics(\%stats);
        die "\n$state\nDESC: $!\nCODE: $stats{ERROR_CODE}";
    }
}

sub ClearError {
    my @stats = qw(
        DATE_MYSQLCREATE     DATE_TRF          DATE_RENUMB      DATE_REDUND
        DATE_PROCLU          DATE_JOINCLUST    DATE_DB_INSERT_REFS
        DATE_DB_INSERT_READS DATE_WRITE_FLANKS DATE_MAP_FLANKS  DATE_MAP_REFFLANKS
        DATE_MAP_INSERT      DATE_EDGES        DATE_INDEX_PCR   DATE_PCR_DUP
        DATE_MAP_DUP         DATE_VNTR_PREDICT DATE_ASSEMBLYREQ DATE_REPORTS
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

    set_statistics({
        ERROR_STEP => 0,
        ERROR_DESC => "",
        ERROR_CODE => 0,
    });
    return 0;
}

sub PrintNextStep {
    # note: this will return 0 regardless if step 0 has been completed or not
    my @stats = qw(
        DATE_MYSQLCREATE     DATE_TRF          DATE_RENUMB      DATE_REDUND
        DATE_PROCLU          DATE_JOINCLUST    DATE_JOINCLUST   DATE_DB_INSERT_REFS
        DATE_DB_INSERT_READS DATE_WRITE_FLANKS DATE_MAP_FLANKS  DATE_MAP_REFFLANKS
        DATE_MAP_INSERT      DATE_EDGES        DATE_INDEX_PCR   DATE_PCR_DUP
        DATE_MAP_DUP         DATE_VNTR_PREDICT DATE_ASSEMBLYREQ DATE_REPORTS
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
            return;
        }
    }
}

sub Stamp {
    print strftime(shift . ": %F %T\n", localtime);
}

sub FinishStep {
    my $name = shift @_;
    my $stats = @_ ? shift @_: {};
    $stats->{"TIME_$name"} = time() - $timestart;

    set_statistics($stats);
    set_datetime("DATE_$name");
    Stamp('End');
    print "\n\n";

    if    ( $STEPEND == $STEP ) { $STEP = 100; }
    elsif ( $CONTINUOUS )       { $STEP++; }
}

####################################

# pipeline error checking

if ( $STEP != 0 ) {

    # clear error?
    if ( $STEP == 100 ) {
        if ( 0 <= $STEPEND && $STEPEND <= 19) { ClearError($STEPEND);}
        else                                  { ClearError();}

        warn "Pipeline error cleared!\n";
        exit;
    }

    my $err_hash = get_statistics(qw(ERROR_STEP ERROR_DESC ERROR_CODE));
    if ( $err_hash->{ERROR_STEP} != 0 ) {
        die "\nPipeline error detected at step $err_hash->{ERROR_STEP} "
            . "(code: $err_hash->{ERROR_CODE}, description: "
            . "'$err_hash->{ERROR_DESC}'). Call this program with step 100 to clear error.\n";
    }

    # get next step
    if ( $STEP == 99 ) {
        PrintNextStep();
        exit;
    }
}

####################################

if ( $STEP == 0 ) {
    print "Executing step #$STEP (creating database).\n";
    Stamp('Start');
    $timestart = time();

    write_sqlite();
    {   no warnings 'once';
        print "SQLite Version: $DBD::SQLite::sqlite_version \n";
    }

    FinishStep('MYSQLCREATE', {
        MAP_ROOT             => $install_dir,
        N_MIN_SUPPORT        => $opts{'MIN_SUPPORT_REQUIRED'},
        MIN_FLANK_REQUIRED   => $opts{'MIN_FLANK_REQUIRED'},
        MAX_FLANK_CONSIDERED => $opts{'MAX_FLANK_CONSIDERED'},
    });
}

if ( $STEP == 1 ) {
    print "Executing step #$STEP (searching for tandem repeats in reads,"
             . " producing profiles and sorting).\n";
    Stamp('Start');
    $timestart = time();

    # Prep
    die "Failed to create output directory $trff.\n"
        unless -r -w -e $trff or mkdir $trff;

    unlink glob("$trff/*.indexhist"), glob("$trff/*.index"),
           glob("$trff/*.leb36"),     glob("$trff/*.reads");

    # Exec
    system("./run_trf_ng.pl",
        $opts{'INPUT_DIR'},
        $trff,
        $config_file,
        $TRF_PARAM,
        $TRF2PROCLU_PARAM,
        $opts{'STRIP_454_KEYTAGS'},
        $opts{'IS_PAIRED_READS'},
        $opts{'NPROCESSES'});
    FlagError('running TRF and TRF2PROCLU');

    FinishStep('TRF', {
        PARAM_TRF       => $TRF_PARAM,
        FOLDER_FASTA    => $opts{'INPUT_DIR'},
        FOLDER_PROFILES => $trff,
    });
}

if ( $STEP == 2 ) {
    print "Executing step #$STEP (reassigning IDs to repeats).\n";
    Stamp('Start');
    $timestart = time();

    unlink glob("$trff/*.leb36.renumbered"), glob("$trff/*.index.renumbered");

    system("./renumber.pl", $trff);
    FlagError('calling renumber.pl on reads profiles folder');

    FinishStep('RENUMB');
}

if ( $STEP == 3 ) {
    print "Executing step #$STEP (eliminating cyclic redundancies).\n";
    Stamp('Start');
    $timestart = time();

    die "Failed to create output directory $processedf.\n"
        unless -r -w -e $processedf or mkdir $processedf;

    unlink glob("$processedf/*.allreads.leb36"), glob("$processedf/*.rotindex");

    # Run redund.exe on read TRs
    system("./redund.exe",
        $trff,
        "$processedf/allreads.leb36",
        "-i");
    FlagError('calling redund.exe on read profiles folder');

    print "Setting additional statistics...\n";
    system("./setdbstats.pl",
        $trff,
        $processedf,
        $config_file);
    FlagError('calling setdbstats.pl');

    FinishStep('REDUND', {
        FOLDER_PROFILES_CLEAN => $processedf,
    });
}

if ( $STEP == 4 ) {
    print "Executing step #$STEP (bipartite clustering of TR profiles).\n";
    Stamp('Start');
    $timestart = time();

    unlink glob("$processedf/*.clu"), glob("$processedf/*.proclu_log");

    system("./checkleb36.pl", $processedf);
    FlagError('checking leb36 files');

    print "Generating temp reference files.\n";
    # 0 for maxerror, means psearch will pick maxerror based on individual flanklength
    # TODO Maybe rewrite run_proclu and psearch to use SQLite db.
    # psearch output can still be files but this step will then collect them and
    # insert into the db. (Maybe merge in next step's `join_clusters.exe`)
    # For now, extract the file we need using the database
    my $dbh = get_ref_dbh( $opts{REFERENCE}, { redo => $opts{REDO_REFDB} } );

    my $reference_folder = File::Temp->newdir();
    # Filtered file, ordered by min representation as given by redund
    # Must negate rids for refset
    my $get_ref_leb36 = q{
        SELECT -rid, length(pattern) AS patsize,
            printf("%.2f", copynum), proflen, proflenrc, profile, profilerc, nA, nC, nG, nT,
            printf("%s|%s", upper(substr(flankleft, -60)), upper(substr(flankright, 0, 61))) AS flanks
        FROM fasta_ref_reps JOIN ref_profiles USING (rid)
            JOIN minreporder USING (rid)
        ORDER BY minreporder.idx ASC};
    my $get_filtered_set_profiles_sth = $dbh->prepare($get_ref_leb36);

    $get_filtered_set_profiles_sth->execute();
    open my $tmp_filt_file_fh, ">", "$reference_folder/reference.leb36";
    while ( my @fields = $get_filtered_set_profiles_sth->fetchrow_array() ) {
        print $tmp_filt_file_fh join( " ", @fields ) . "\n";
    }
    close $tmp_filt_file_fh;

    # Get rotindex saved in db
    my ($rotindex_str) = $dbh->selectrow_array("SELECT rotindex FROM files");
    die "Error getting rotindex file from filtered set. Try rerunning with --redo option\n"
        unless ($rotindex_str);
    $dbh->disconnect();

    # Need to negate all indices
    open my $tmp_rotindex, ">", "$reference_folder/reference.leb36.rotindex";
    $rotindex_str =~ s/(\d+)/-$1/g;
    print $tmp_rotindex $rotindex_str;
    close $tmp_rotindex;

    print "Running proclu.\n";
    system("./run_proclu.pl",
        1,
        $processedf,
        $reference_folder,
        $CLUST_PARAMS,
        $opts{'NPROCESSES'},
        $PROCLU_EXECUTABLE,
        0,
        $opts{'MAX_FLANK_CONSIDERED'});
    FlagError('performing bipartite clustering of tandem repeats profiles failed');

    # remove clusters with no references
    print "Removing clusters with no references (clara/pam split).\n";
    opendir( DIR, $processedf );
    my @files = grep( /clu$/, readdir(DIR) ); # KA: more file grepping
    closedir(DIR);

    foreach my $file (@files) {
        if ( open( my $cluf, "<", "$processedf/$file" ) ) {
            open( my $clufout, ">",
                "$processedf/$file.clean" );
            while (<$cluf>) {
                if (/-/) { print $clufout $_; }
            }
            close($cluf);
            close($clufout);
            rename "$processedf/$file.clean", "$processedf/$file";
        }
    }

    FinishStep('PROCLU', {
        FILE_REFERENCE_LEB => $opts{'REFERENCE'} . ".leb36",
        FILE_REFERENCE_SEQ => $opts{'REFERENCE'} . ".seq",
        PARAM_PROCLU => $CLUST_PARAMS,
    });
}

if ( $STEP == 5 ) {
    print "Executing step #$STEP (joining clusters from different proclu runs on reference ids).\n";
    Stamp('Start');
    $timestart = time();

    unlink "$processedf/allwithdups.clusters";

    system("./join_clusters.exe",
        $processedf,
        "$processedf/allwithdups.clusters");
    FlagError('joining clusters from different proclu runs on reference ids');

    FinishStep('JOINCLUST');
}

# Steps 6 and 7 merged or removed
if ( $STEPEND == 6 or $STEPEND == 7 ) { $STEP = 100; }
elsif ( $STEP < 8 and $CONTINUOUS )   { $STEP = 8; }

if ( $STEP == 8 ) {
    print "Executing step #$STEP (inserting READS flanks into database).\n";
    Stamp('Start');
    $timestart = time();

    system("./insert_reads.pl",
        "$processedf/allwithdups.clusters",
        $trff,
        $opts{'STRIP_454_KEYTAGS'},
        $config_file);
    FlagError('inserting READS flanks into database');

    FinishStep('DB_INSERT_READS');
}

if ( $STEP == 9 ) {
    print "Executing step #$STEP (outputting flanks inside each cluster).\n";
    Stamp('Start');
    $timestart = time();

    unlink "$processedf/allwithdups.flanks";

    system("./run_flankcomp.pl",
        "$processedf/allwithdups.clusters",
        $config_file,
        File::Temp->newdir(),
        "$processedf/allwithdups.flanks");
    FlagError('outputting flanks inside each cluster');

    FinishStep('WRITE_FLANKS');
}

if ( $STEP == 10 ) {
    print "Executing step #$STEP (aligning ref-read flanks).\n";
    Stamp('Start');
    $timestart = time();

    my $outf = "$processedf/out";
    die "Failed to create output directory $outf.\n"
        unless -r -w -e $outf or mkdir $outf;
    unlink glob "$outf/*";

    # 0 for maxerror, means flankalign will pick maxerror based on individual flanklength
    system("./flankalign.exe",
        $outf,
        "$processedf/result",
        "$processedf/allwithdups.flanks",
        0,
        $opts{'MAX_FLANK_CONSIDERED'},
        $opts{'NPROCESSES'},
        15);
    FlagError('aligning ref-read flanks');

    FinishStep('MAP_FLANKS');
}

# Step 11 merged or removed
if ( $STEPEND == 11 )                { $STEP = 100; }
elsif ( $STEP < 12 and $CONTINUOUS ) { $STEP = 12; }

if ( $STEP == 12 ) {
    print "Executing step #$STEP (inserting map and rankflank information into database).\n";
    Stamp('Start');
    $timestart = time();

    system("./run_rankflankmap.pl",
        "$processedf/allwithdups.clusters",
        "$processedf/out",
        $config_file);
    FlagError('inserting map and rankflank information into database');

    FinishStep('MAP_INSERT');
}

if ( $STEP == 13 ) {
    print "Executing step #$STEP (calculating edges).\n";
    Stamp('Start');
    $timestart = time();

    system("./run_edges.pl",
        $opts{'REFERENCE'},
        File::Temp->newdir(), # edges folder
        $config_file,
        $MINPROFSCORE,
        $opts{'NPROCESSES'},
        $PROCLU_EXECUTABLE);
    FlagError('calculating edges');

    FinishStep('EDGES');
}

if ( $STEP == 14 ) {
    print "Executing step #$STEP (calculating PCR duplicates).\n";
    Stamp('Start');
    $timestart = time();

    if ( $opts{KEEPPCRDUPS} ) {
        print "Notice: not deleting duplicates.\n";
    }

    my $bestf = "$processedf/best";
    die "Failed to create output directory $bestf.\n"
        unless -r -w -e $bestf or mkdir $bestf;
    unlink glob "$bestf/*";
    #unlink "$output_folder/$opts{RUN_NAME}.pcr_dup.txt";
    #unlink "$output_folder/$opts{RUN_NAME}.ties.txt";
    #unlink "$output_folder/$opts{RUN_NAME}.ties_entries.txt";

    system("./pcr_dup.pl",
        "$processedf/best",
        $output_folder,
        $opts{'RUN_NAME'},
        $config_file,
        $opts{'NPROCESSES'},
        $opts{'KEEPPCRDUPS'});
    FlagError('calculating PCR duplicates');

    remove_tree("$processedf/best", {safe => 1});

    FinishStep('PCR_DUP');
}

# Step 15 merged into 14
if ( $STEPEND == 15 )                { $STEP = 100; }
elsif ( $STEP < 16 and $CONTINUOUS ) { $STEP = 16; }

if ( $STEP == 16 ) {
    print "Executing step #$STEP (removing mapped duplicates).\n";
    Stamp('Start');
    $timestart = time();

    #unlink "$output_folder/$opts{RUN_NAME}.map_dup.txt";

    system("./map_dup.pl",
        $config_file,
        "$output_folder/$opts{RUN_NAME}.map_dup.txt"); # file writing turned off
    FlagError('calculating mapped duplicates failed');

    FinishStep('MAP_DUP');
}

if ( $STEP == 17 ) {
    print "Executing step #$STEP (computing variability).\n";
    Stamp('Start');
    $timestart = time();

    system("./run_variability.pl",
        "$processedf/allwithdups.clusters",
        $config_file,
        $opts{'MIN_FLANK_REQUIRED'});
    FlagError('computing variability');

    FinishStep('VNTR_PREDICT');
}

# Step 18 merged or removed
if ( $STEPEND == 18 )                { $STEP = 100; }
elsif ( $STEP < 19 and $CONTINUOUS ) { $STEP = 19; }

if ( $STEP == 19 ) {
    print "Executing step #$STEP (final database update).\n";
    Stamp('Start');
    $timestart = time();

    unlink glob "$output_folder/*.vcf";

    # further output files (vcfs)
    system("./updaterefs.pl",
        $opts{'RUN_NAME'},
        $config_file,
        "$output_folder/$opts{RUN_NAME}",
        $VERSION);
    FlagError('final database update');

    FinishStep('REPORTS');
}

if ( $STEP == 20 ) {
    print "Executing step #$STEP (database minify and cleanup).\n";
    Stamp('Start');
    $timestart = time();

    # Reduce Database
    system('./reduce_db.pl',
        $opts{'RUN_NAME'},
        $output_folder,
        $opts{'REFERENCE'} . ".db",
        $opts{'READ_LENGTH'},
        File::Temp->newdir(),
        $opts{'NPROCESSES'});
    FlagError('db conversion and compression');

    # Cleanup
    print "File Cleanup Time!\n";
    remove_tree("$processedf/out", {safe => 1});
    if ($opts{'CLEANUP'}) {
        remove_tree($trff, {safe => 1});
        remove_tree($processedf, {safe => 1});
    }

    my $sn = 'DB_CONVERSION_AND_READ_COMPRESSION';
    FinishStep($sn);

    # Migrate final stats to reduced database
    my $dbf  = "$output_folder/$opts{RUN_NAME}.db";
    my $dbf2 = "$output_folder/$opts{RUN_NAME}_rl$opts{READ_LENGTH}.db";
    my $dbh = DBI->connect("DBI:SQLite:dbname=$dbf");
    $dbh->do(qq(ATTACH DATABASE "$dbf2" as newdb));
    $dbh->do(q{
    UPDATE newdb.stats
    SET TIME_$sn = (select TIME_$sn from stats),
        DATE_$sn = (select DATE_$sn from stats)});
}

print "Finished!\n";
