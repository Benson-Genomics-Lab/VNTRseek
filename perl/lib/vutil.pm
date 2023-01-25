package vutil;
use v5.24;
use Cwd;
use Cwd 'abs_path';
use DBI;
use Carp;
use FindBin;
use File::Temp;
use Try::Tiny;
use POSIX 'strftime';


if ( $ENV{DEBUG} ) {
    use Data::Dumper;
}

use base 'Exporter';
our @EXPORT_OK = qw(trim read_config_file get_config validate_config
    get_dbh get_ref_dbh make_refseq_db load_refprofiles_db
    run_redund write_sqlite set_statistics get_statistics set_datetime
    print_config trim sqlite_install_RC_function
    gen_exec_array_cb vs_db_insert);

# vutil.pm
# author: Yevgeniy Gelfand, Yozen Hernandez
# maintainer: Kyler Anderson
# create date: Oct 30, 2010
# last edited: Jul 29, 2022
# function: common functions for vntr pipleline scripts.
# Provide functions for database management, configuration
# reading/writing, and miscellaneous common utilities.

my %cnf = (); # Run Configuation
my %valid_stats = (
    MAP_ROOT                                            => 1,
    PARAM_TRF                                           => 1,
    PARAM_PROCLU                                        => 1,
    FOLDER_FASTA                                        => 1,
    FOLDER_PROFILES                                     => 1,
    FOLDER_PROFILES_CLEAN                               => 1,
    FOLDER_REFERENCE                                    => 1,
    FILE_REFERENCE_LEB                                  => 1,
    FILE_REFERENCE_SEQ                                  => 1,
    NUMBER_READS                                        => 1,
    NUMBER_TRS_IN_READS                                 => 1,
    NUMBER_TRS_IN_READS_GE7                             => 1,
    NUMBER_READS_WITHTRS                                => 1,
    NUMBER_READS_WITHTRS_GE7                            => 1,
    NUMBER_READS_WITHTRS_GE7_AFTER_REDUND               => 1,
    NUMBER_TRS_IN_READS_AFTER_REDUND                    => 1,
    NUMBER_REF_TRS                                      => 1,
    NUMBER_REFS_TRS_AFTER_REDUND                        => 1,
    CLUST_NUMBER_OF_PROCLU_CLUSTERS                     => 1,
    CLUST_NUMBER_OF_PROCLU_CLUSTERS_BEFORE_REJOIN       => 1,
    CLUST_NUMBER_OF_EXACTPAT_CLUSTERS                   => 1,
    CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS                => 1,
    CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS               => 1,
    CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER       => 1,
    CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER      => 1,
    CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER            => 1,
    CLUST_LARGEST_NUMBER_OF_TRS_IN_EXACTPAT_CLUSTER     => 1,
    CLUST_LARGEST_NUMBER_OF_REFS_IN_EXACTPAT_CLUSTER    => 1,
    CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR            => 1,
    CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR        => 1,
    NUMBER_REFS_VNTR_SPAN_N                             => 1,
    NUMBER_REFS_SINGLE_REF_CLUSTER                      => 1,
    NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED    => 1,
    NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED => 1,
    NUMBER_MAPPED                                       => 1,
    NUMBER_RANK                                         => 1,
    NUMBER_RANKFLANK                                    => 1,
    INTERSECT_RANK_AND_RANKFLANK                        => 1,
    INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR             => 1,
    BBB_WITH_MAP_DUPS                                   => 1,
    BBB                                                 => 1,
    RANK_EDGES_OVERCUTOFF                               => 1,
    RANK_REMOVED_SAMEREF                                => 1,
    RANK_REMOVED_SAMESEQ                                => 1,
    RANK_REMOVED_PCRDUP                                 => 1,
    RANKFLANK_EDGES_INSERTED                            => 1,
    RANKFLANK_REMOVED_SAMEREF                           => 1,
    RANKFLANK_REMOVED_SAMESEQ                           => 1,
    RANKFLANK_REMOVED_PCRDUP                            => 1,
    TIME_MYSQLCREATE                                    => 1,
    TIME_TRF                                            => 1,
    TIME_RENUMB                                         => 1,
    TIME_REDUND                                         => 1,
    TIME_PROCLU                                         => 1,
    TIME_JOINCLUST                                      => 1,
    TIME_DB_INSERT_REFS                                 => 1,
    TIME_DB_INSERT_READS                                => 1,
    TIME_WRITE_FLANKS                                   => 1,
    TIME_MAP_FLANKS                                     => 1,
    TIME_MAP_REFFLANKS                                  => 1,
    TIME_MAP_INSERT                                     => 1,
    TIME_EDGES                                          => 1,
    TIME_INDEX_PCR                                      => 1,
    TIME_PCR_DUP                                        => 1,
    TIME_MAP_DUP                                        => 1,
    TIME_VNTR_PREDICT                                   => 1,
    TIME_ASSEMBLYREQ                                    => 1,
    TIME_REPORTS                                        => 1,
    TIME_DB_CONVERSION_AND_READ_COMPRESSION             => 1,
    DATE_MYSQLCREATE                                    => 1,
    DATE_TRF                                            => 1,
    DATE_RENUMB                                         => 1,
    DATE_REDUND                                         => 1,
    DATE_PROCLU                                         => 1,
    DATE_JOINCLUST                                      => 1,
    DATE_DB_INSERT_REFS                                 => 1,
    DATE_DB_INSERT_READS                                => 1,
    DATE_WRITE_FLANKS                                   => 1,
    DATE_MAP_FLANKS                                     => 1,
    DATE_MAP_REFFLANKS                                  => 1,
    DATE_MAP_INSERT                                     => 1,
    DATE_EDGES                                          => 1,
    DATE_INDEX_PCR                                      => 1,
    DATE_PCR_DUP                                        => 1,
    DATE_MAP_DUP                                        => 1,
    DATE_VNTR_PREDICT                                   => 1,
    DATE_ASSEMBLYREQ                                    => 1,
    DATE_REPORTS                                        => 1,
    DATE_DB_CONVERSION_AND_READ_COMPRESSION             => 1,
    ERROR_STEP                                          => 1,
    ERROR_DESC                                          => 1,
    ERROR_CODE                                          => 1,
    N_MIN_SUPPORT                                       => 1,
    MIN_FLANK_REQUIRED                                  => 1,
    MAX_FLANK_CONSIDERED                                => 1,
);

################################################################
=item C<trim>

Takes an input string and removes the whitespace from the start
and end of a string.

=cut

sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

################################################################
sub read_config_file {
    croak "read_config_file() expects 1 parameter.\n"
        unless (@_);

    my $file_loc = shift;
    my %opts = ();

    # read config
    open(my $cnfh, "<", "$file_loc") or die "$!";
    while (<$cnfh>) {
        chomp;
        if (/^\#/) { next;} # skip comments

        if (/(.+)=(.*)/) {
            my $key = trim($1);
            my $val = trim($2);
            $val =~ s/\s*\#.*//;    # strip end comments
            $opts{ uc($key) } = $val;
        }
    }
    close($cnfh);
    return %opts;
}

################################################################
sub get_config {
    croak "get_config() expects an even number of parameters\n"
        unless (scalar @_ % 2 == 0);
    my %args = @_;

    # load default settings
    my %opts = read_config_file("$FindBin::RealBin/defaults.vs.cnf");

    # load local settings
    if (exists $args{'CONFIG'} and $args{'CONFIG'}) {
        croak "Cannot access local config at $args{'CONFIG'}. Please check name and permissions.\n"
            unless -r -e $args{'CONFIG'};

        my %local = read_config_file($args{'CONFIG'});
        for my $key (keys %local) {
            $opts{$key} = $local{$key};
        }
    }

    # write command line settings
    for my $key (keys %args) {
        $opts{$key} = $args{$key};
    }

    # Backwards compatibility option substitutions
    my %renames = (
        'FASTA_DIR', 'INPUT_DIR',
        'OUTPUT_ROOT', 'OUTPUT_DIR',
        'DBSUFFIX', 'RUN_NAME'
    );

    for my $old_name (keys %renames) {
        if ($opts{$old_name}) { $opts{$renames{$old_name}} = $opts{$old_name};}
    }

    # Since all of the sub-scripts are called in a subprocess,
    #   we consume a process before ever forking.
    # To avoid using more processors than advertised,
    #   we need to shave one off.
    $opts{'NPROCESSES'} -= 1;

    # Resolve relative paths for future chdirs
    for my $file ('INPUT_DIR', 'OUTPUT_DIR') {
        if ($opts{$file}) { $opts{$file} = abs_path($opts{$file});}
    }

    # keep a local pointer and send back the results
    %cnf = %opts;
    return %opts;
}

################################################################
sub validate_config {
    my $please_check = "Please check the command line arguments or the configuration file.\n";

    # Validation
    croak("Please set the run name (RUN_NAME) on the command line or in the configuration file.\n")
        unless $cnf{'RUN_NAME'};

    croak("Please set the input directory (INPUT_DIR) on the command line or in the configuration file.\n")
        unless $cnf{'INPUT_DIR'};
    croak("Input directory '$cnf{'INPUT_DIR'}' could not be read. Check the spelling and/or your permissions.")
        unless -r -e $cnf{'INPUT_DIR'};

    croak("OUTPUT_DIR value cannot be blank. " . $please_check)
        unless $cnf{'OUTPUT_DIR'};
    croak("Output directory '$cnf{'OUTPUT_DIR'}' could not be accessed. Check the spelling and/or your permissions.")
        unless (-e $cnf{'OUTPUT_DIR'} or mkdir( $cnf{'OUTPUT_DIR'} )) and -x -w -r $cnf{'OUTPUT_DIR'};

    croak("SERVER value cannot be blank. " . $please_check)
        unless $cnf{'SERVER'};


    croak("PLOIDY value must be > 0. " . $please_check)
        unless $cnf{'PLOIDY'} > 0;

    croak("READ_LENGTH value must be > 0. " . $please_check)
        unless $cnf{'READ_LENGTH'} > 0;

    croak("IS_PAIRED_READS value must be 0/1. " . $please_check)
        unless $cnf{'IS_PAIRED_READS'} >= 0;

    croak("STRIP_454_KEYTAGS value must be 0/1. " . $please_check)
        unless $cnf{'STRIP_454_KEYTAGS'} >= 0;

    croak("KEEPPCRDUPS value must be 0/1. " . $please_check)
        unless $cnf{'KEEPPCRDUPS'} >= 0;


    croak("MIN_FLANK_REQUIRED value must be > 0. " . $please_check)
        unless $cnf{'MIN_FLANK_REQUIRED'} > 0;

    croak("MAX_FLANK_CONSIDERED value must be > 0. " . $please_check)
        unless $cnf{'MAX_FLANK_CONSIDERED'} > 0;

    croak("MIN_SUPPORT_REQUIRED value must be > 0. " . $please_check)
        unless $cnf{'MIN_SUPPORT_REQUIRED'} > 0;


    croak("NPROCESSES value must be > 1. " . $please_check)
        unless $cnf{'NPROCESSES'} > 0; # Intentional discrepancy, see get_config


    croak("Please set the reference base name (REFERENCE) on the command line or in the configuration file.\n")
        unless $cnf{'REFERENCE'};

    croak("REFERENCE_INDIST_PRODUCE value must be 0/1. " . $please_check)
        unless $cnf{'REFERENCE_INDIST_PRODUCE'} >= 0;

    if ( !( -e $cnf{'REFERENCE'} . ".db" ) || $cnf{'REDO_REFDB'} ) {
        for my $ext ('.leb36', '.seq', '.indist') {
            croak("Reference file '$cnf{'REFERENCE'}$ext' not found!\n"
                  . "You may need to rebuild the reference db with REDO_REFDB option.\n")
                unless -r -e $cnf{'REFERENCE'} . $ext;
        }
    }
}


################################################################
sub print_config {
    my $output_folder = "$cnf{'OUTPUT_DIR'}/vntr_$cnf{'RUN_NAME'}";
    my $config_file = "$output_folder/$cnf{'RUN_NAME'}.vs.cnf";
    my $when = strftime("%c", localtime);
    my $padproc = $cnf{'NPROCESSES'} + 1;
    open( my $config_fh, ">", $config_file )
        or die "Cannot open '$config_file' for writing!\n";

    print $config_fh <<CNF;
# VNTRseek complete configuration
# $when
# You can provide this file with the --CONFIG
#   option to reproduce the latest run.
# This will always be generated here and
#   overwritten to match the last configuration
#   with this output location.
# Use vntrseek --HELP for option descriptions.

RUN_NAME=$cnf{"RUN_NAME"}
INPUT_DIR=$cnf{"INPUT_DIR"}
OUTPUT_DIR=$cnf{"OUTPUT_DIR"}
SERVER=$cnf{"SERVER"}

REFERENCE=$cnf{"REFERENCE"}
REDO_REFDB=$cnf{"REDO_REFDB"}
REFERENCE_INDIST_PRODUCE=$cnf{"REFERENCE_INDIST_PRODUCE"}

PLOIDY=$cnf{"PLOIDY"}
IS_PAIRED_READS=$cnf{"IS_PAIRED_READS"}
READ_LENGTH=$cnf{"READ_LENGTH"}
KEEPPCRDUPS=$cnf{"KEEPPCRDUPS"}
STRIP_454_KEYTAGS=$cnf{"STRIP_454_KEYTAGS"}

MIN_FLANK_REQUIRED=$cnf{"MIN_FLANK_REQUIRED"}
MAX_FLANK_CONSIDERED=$cnf{"MAX_FLANK_CONSIDERED"}
MIN_SUPPORT_REQUIRED=$cnf{"MIN_SUPPORT_REQUIRED"}

NPROCESSES=$padproc
CNF

    close($config_fh);
    chmod 0640, $config_file; #why?
    return $config_file;
}

####################################
sub set_statistics {

# my $argc = @_;
# if ( $argc < 2 ) {
#     croak "set_statistics: expects at least 2 parameters, passed $argc !\n";
# }

    my $stats = shift;
    my $dbh   = get_dbh();
    my ( $sql_clause, @sql_qual, @sql_bind );
    ( $ENV{DEBUG} ) && warn Dumper($stats) . "\n";

    while ( my ( $key, $val ) = each %$stats ) {
        next unless exists $valid_stats{$key};
        ( $ENV{DEBUG} ) && warn "Setting stat: $key to $val\n";
        push @sql_qual, "$key=?";
        push @sql_bind, $val;
    }

    $sql_clause = join ",", @sql_qual;
    my $sth = $dbh->prepare("UPDATE stats SET $sql_clause")
        or croak "Couldn't prepare statement: " . $dbh->errstr;

    $sth->execute(@sql_bind)    # Execute the query
        or croak "Couldn't execute statement: " . $sth->errstr;

    $dbh->disconnect();
}

################################################################

sub set_datetime {

    my $argc = @_;
    if ( $argc < 1 ) {
        croak "set_datetime: expects 1 parameter, passed $argc !\n";
    }

    my $NAME  = shift;
    my $VALUE = strftime( "%F %T", localtime );

    return set_statistics( { $NAME, $VALUE } );
}

####################################

##
## @brief      Gets the statistics.
##
## @return     The statistics.
##
sub get_statistics {

  # my $argc = @_;
  # if ( $argc < 2 ) {
  #     die "get_statistics: expects at least 2 parameters, passed $argc !\n";
  # }

    # my $RUN_NAME = shift;
    my @stats = grep { exists $valid_stats{$_} } @_;
    my $dbh = get_dbh( { readonly => 1 } );
    my $sql_clause;
    ( $ENV{DEBUG} ) && warn Dumper( \@stats ) . "\n";

    $sql_clause = join ", ", @stats;
    ( $ENV{DEBUG} ) && warn "Getting stats: " . $sql_clause . "\n";

    my $sql_res = $dbh->selectrow_hashref("SELECT $sql_clause FROM stats");

    $dbh->disconnect();
    return $sql_res;
}

################################################################
sub _init_ref_dbh {
    my ( $refbase, $opts ) = @_;
    my $refdbfile = $refbase . ".db";
    my $refseq    = $refbase . ".seq";
    my $refleb36  = $refbase . ".leb36";
    my $refindist = $refbase . ".indist";
    my %dbi_opts  = (
        AutoCommit                 => 1,
        RaiseError                 => 1,
        PrintError                 => 0,
        sqlite_see_if_its_a_number => 1,
    );

    # TODO Maybe first connect to a temp location, backup database
    # to that location then return handle to that location. This
    # is primarily for running on clusters/over NFS.
    my $dbh = DBI->connect( "DBI:SQLite:dbname=" . $refdbfile,
        undef, undef, \%dbi_opts )
        or die "Could not connect to database "
        . $refdbfile
        . ": $DBI::errstr";

    # Set some pragmas to make this part faster
    $dbh->do(q{PRAGMA synchronous = OFF});

    # $dbh->do(q{PRAGMA journal_mode = WAL});

    unless ( exists $opts->{skip_refseq} ) {
        make_refseq_db( $dbh, $refseq, $opts->{redo} );
    }

    load_refprofiles_db( $dbh, $refleb36, $opts->{redo} );

    if ( -e $refindist ) {
        set_indist( $dbh, $refindist, $opts->{redo} );
    }

    if ( $opts->{redo} ) {
        $cnf{REDO_REFDB} = 0;
    }

    # END TODO
    $dbh->disconnect();

    # Return the db handle
    return DBI->connect(
        "DBI:SQLite:dbname=" . $refdbfile,
        undef, undef,
        {   %dbi_opts,
            ReadOnly => ( exists $opts->{readonly} ) ? $opts->{readonly} : 0
        }
    );
}
################################################################
sub get_ref_dbh {
    croak "Error: expecting at least 1 argument, got none.\n"
        unless ( @_ >= 1 );

    my ( $refbase, $opts ) = @_;
    if ( $ENV{DEBUG} ) { warn "get_ref_dbh arguments: " . Dumper( \@_ ) . "\n" };
    my $dbh = _init_ref_dbh( $refbase, $opts );

    return $dbh;
}

################################################################
sub get_dbh {

    # Connection options
    my $opts = shift // {};

    my $dbh;
    my $dbfile
        = "$cnf{OUTPUT_DIR}/vntr_$cnf{RUN_NAME}/$cnf{RUN_NAME}.db";

    # If the db file doesn't exist, or is size 0, initialize the db
    if ( !-e $dbfile || -z $dbfile ) {
        write_sqlite();
    }

    my %dbi_opts = (
        AutoCommit                 => 1,
        RaiseError                 => 1,
        PrintError                 => 0,
        sqlite_see_if_its_a_number => 1,
        ReadOnly => 1 * ( exists $opts->{readonly} && $opts->{readonly} )
    );

    $dbh = DBI->connect( "DBI:SQLite:dbname=$dbfile", undef, undef, \%dbi_opts )
        or die "Could not connect to database $dbfile: $DBI::errstr";

    my ($stats_schema) = $dbh->selectrow_array(
        q{SELECT sql FROM sqlite_master WHERE name = 'stats'}
    );

    # If the stats table does not exist, init the db.
    unless ($stats_schema) {
        write_sqlite();
    }

    my ($schema_ver) = $dbh->selectrow_array(q{PRAGMA user_version});

    if ( $schema_ver < 1 ) {
        $dbh->disconnect();
        $dbh = DBI->connect( "DBI:SQLite:dbname=$dbfile", undef, undef,
            { %dbi_opts, ReadOnly => 0 } )
            or die "Could not connect to database $dbfile: $DBI::errstr";
        $dbh->do(q{PRAGMA foreign_keys = off});
        $dbh->begin_work();
        $dbh->do(q{ALTER TABLE stats RENAME TO _old_stats});
        $dbh->do(
            q{CREATE TABLE `stats` (
  `id` integer NOT NULL,
  `MAP_ROOT` varchar(500) DEFAULT NULL,
  `PARAM_TRF` varchar(500) DEFAULT NULL,
  `PARAM_PROCLU` varchar(500) DEFAULT NULL,
  `FOLDER_FASTA` varchar(500) DEFAULT NULL,
  `FOLDER_PROFILES` varchar(500) DEFAULT NULL,
  `FOLDER_PROFILES_CLEAN` varchar(500) DEFAULT NULL,
  `FOLDER_REFERENCE` varchar(500) DEFAULT NULL,
  `FILE_REFERENCE_LEB` varchar(500) DEFAULT NULL,
  `FILE_REFERENCE_SEQ` varchar(500) DEFAULT NULL,
  `NUMBER_READS` integer DEFAULT NULL,
  `NUMBER_TRS_IN_READS` integer DEFAULT NULL,
  `NUMBER_TRS_IN_READS_GE7` integer DEFAULT NULL,
  `NUMBER_READS_WITHTRS` integer DEFAULT NULL,
  `NUMBER_READS_WITHTRS_GE7` integer DEFAULT NULL,
  `NUMBER_READS_WITHTRS_GE7_AFTER_REDUND` integer DEFAULT NULL,
  `NUMBER_TRS_IN_READS_AFTER_REDUND` integer DEFAULT NULL,
  `NUMBER_REF_TRS` integer DEFAULT NULL,
  `NUMBER_REFS_TRS_AFTER_REDUND` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_PROCLU_CLUSTERS` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_PROCLU_CLUSTERS_BEFORE_REJOIN` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_EXACTPAT_CLUSTERS` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_TRS_IN_EXACTPAT_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_REFS_IN_EXACTPAT_CLUSTER` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR` integer DEFAULT NULL,
  `NUMBER_REFS_VNTR_SPAN_N` integer DEFAULT NULL,
  `NUMBER_REFS_SINGLE_REF_CLUSTER` integer DEFAULT NULL,
  `NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED` integer DEFAULT NULL,
  `NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED` integer DEFAULT NULL,
  `NUMBER_MAPPED` integer DEFAULT NULL,
  `NUMBER_RANK` integer DEFAULT NULL,
  `NUMBER_RANKFLANK` integer DEFAULT NULL,
  `INTERSECT_RANK_AND_RANKFLANK` integer DEFAULT NULL,
  `INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR` integer DEFAULT NULL,
  `BBB_WITH_MAP_DUPS` integer DEFAULT NULL,
  `BBB` integer DEFAULT NULL,
  `RANK_EDGES_OVERCUTOFF` integer DEFAULT NULL,
  `RANK_REMOVED_SAMEREF` integer DEFAULT NULL,
  `RANK_REMOVED_SAMESEQ` integer DEFAULT NULL,
  `RANK_REMOVED_PCRDUP` integer DEFAULT NULL,
  `RANKFLANK_EDGES_INSERTED` integer DEFAULT NULL,
  `RANKFLANK_REMOVED_SAMEREF` integer DEFAULT NULL,
  `RANKFLANK_REMOVED_SAMESEQ` integer DEFAULT NULL,
  `RANKFLANK_REMOVED_PCRDUP` integer DEFAULT NULL,
  `TIME_MYSQLCREATE` integer DEFAULT NULL,
  `TIME_TRF` integer DEFAULT NULL,
  `TIME_RENUMB` integer DEFAULT NULL,
  `TIME_REDUND` integer DEFAULT NULL,
  `TIME_PROCLU` integer DEFAULT NULL,
  `TIME_JOINCLUST` integer DEFAULT NULL,
  `TIME_DB_INSERT_REFS` integer DEFAULT NULL,
  `TIME_DB_INSERT_READS` integer DEFAULT NULL,
  `TIME_WRITE_FLANKS` integer DEFAULT NULL,
  `TIME_MAP_FLANKS` integer DEFAULT NULL,
  `TIME_MAP_REFFLANKS` integer DEFAULT NULL,
  `TIME_MAP_INSERT` integer DEFAULT NULL,
  `TIME_EDGES` integer DEFAULT NULL,
  `TIME_INDEX_PCR` integer DEFAULT NULL,
  `TIME_PCR_DUP` integer DEFAULT NULL,
  `TIME_MAP_DUP` integer DEFAULT NULL,
  `TIME_VNTR_PREDICT` integer DEFAULT NULL,
  `TIME_ASSEMBLYREQ` integer DEFAULT NULL,
  `TIME_REPORTS` integer DEFAULT NULL,
  `TIME_DB_CONVERSION_AND_READ_COMPRESSION` integer DEFAULT NULL,
  `DATE_MYSQLCREATE` text DEFAULT NULL,
  `DATE_TRF` text DEFAULT NULL,
  `DATE_RENUMB` text DEFAULT NULL,
  `DATE_REDUND` text DEFAULT NULL,
  `DATE_PROCLU` text DEFAULT NULL,
  `DATE_JOINCLUST` text DEFAULT NULL,
  `DATE_DB_INSERT_REFS` text DEFAULT NULL,
  `DATE_DB_INSERT_READS` text DEFAULT NULL,
  `DATE_WRITE_FLANKS` text DEFAULT NULL,
  `DATE_MAP_FLANKS` text DEFAULT NULL,
  `DATE_MAP_REFFLANKS` text DEFAULT NULL,
  `DATE_MAP_INSERT` text DEFAULT NULL,
  `DATE_EDGES` text DEFAULT NULL,
  `DATE_INDEX_PCR` text DEFAULT NULL,
  `DATE_PCR_DUP` text DEFAULT NULL,
  `DATE_MAP_DUP` text DEFAULT NULL,
  `DATE_VNTR_PREDICT` text DEFAULT NULL,
  `DATE_ASSEMBLYREQ` text DEFAULT NULL,
  `DATE_REPORTS` text DEFAULT NULL,
  `DATE_DB_CONVERSION_AND_READ_COMPRESSION` text DEFAULT NULL,
  `ERROR_STEP` integer NOT NULL DEFAULT '0',
  `ERROR_DESC` varchar(500) NOT NULL DEFAULT '',
  `ERROR_CODE` integer NOT NULL DEFAULT '0',
  `N_MIN_SUPPORT` integer NOT NULL DEFAULT '0',
  `MIN_FLANK_REQUIRED` integer NOT NULL DEFAULT '0',
  `MAX_FLANK_CONSIDERED` integer NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`)
);}
        );
        $dbh->do(q{INSERT INTO stats SELECT * FROM _old_stats});
        $dbh->do(q{DROP TABLE _old_stats});
        $dbh->commit();
        $dbh->do(q{PRAGMA foreign_keys = on});
        $dbh->do(q{PRAGMA user_version = 1});
        $schema_ver = 1;
        $dbh->disconnect();
        $dbh = DBI->connect( "DBI:SQLite:dbname=$dbfile",
            undef, undef, \%dbi_opts )
            or die "Could not connect to database $dbfile: $DBI::errstr";
    }

    # Set default journal to write-ahead log
    # $dbh->do(q{PRAGMA journal_mode = WAL})
    # or die "Could not do statement on database $dbfile: $DBI::errstr";

    # 800MB cache
    $dbh->do("PRAGMA cache_size = 800000")
        or die "Could not do statement on database $dbfile: $DBI::errstr";

    # Attach reference set database
    if ( exists $opts->{userefdb} && $opts->{userefdb} ) {
        my $refdbfile = $cnf{REFERENCE} . ".db";

        # First initialize the reference DB, if needed
        # Don't save the returned dbh
        # _init_ref_dbh( $cnf{REFERENCE},
        #     { redo => $cnf{REDO_REFDB} } );
        ( $ENV{DEBUG} )
            && warn "Connecting to refseq db at $refdbfile\n";
        try {
            $dbh->do(qq{ATTACH DATABASE "$refdbfile" AS refdb});
        }
        catch {
            warn "$_\n";
            if (/unable to open database/) {
                die "Try running VNTRseek as:\n"
                    . "\n\tvntrseek --reference $cnf{REFERENCE} --redo_refdb\n\n"
                    . "then clear the error and try again.\n";
            }
            else {
                exit 1;
            }
        }
    }

    return $dbh;
}

####################################
sub make_refseq_db {

    # TODO DBI error checking
    unless ( @_ == 3 ) {
        croak "Error: expecting 3 arguments, got " . @_ * 1 . "\n";
    }
    my ( $dbh, $refseq, $redo ) = @_;

    ( $ENV{DEBUG} ) && warn "CWD = " . getcwd() . "\n";

    # Create the table in this new db.
    my $create_fasta_ref_reps_q
        = q{CREATE TABLE IF NOT EXISTS `fasta_ref_reps` (
    `rid` integer NOT NULL,
    `firstindex` integer NOT NULL,
    `lastindex` integer NOT NULL,
    `copynum` float NOT NULL,
    `head` varchar(500) DEFAULT NULL,
    `flankleft` text COLLATE BINARY,
    `pattern` text NOT NULL,
    `sequence` text NOT NULL,
    `flankright` text COLLATE BINARY,
    `conserved` float DEFAULT NULL,
    `comment` varchar(500) DEFAULT NULL,
    `is_singleton` integer NOT NULL DEFAULT '1',
    `is_dist` integer NOT NULL DEFAULT '1',
    `is_indist` integer NOT NULL DEFAULT '0',
    PRIMARY KEY (`rid`),
    UNIQUE (`rid`,`comment`))};
    my $fasta_ref_reps_index_q = q{CREATE INDEX IF NOT EXISTS
        "idx_fasta_ref_reps_head" ON "fasta_ref_reps" (`head`)};
    my $fasta_ref_reps_patsize_q = q{CREATE INDEX IF NOT EXISTS
        "idx_fasta_ref_reps_patsize" ON
        "fasta_ref_reps" (LENGTH(pattern))};
    my $fasta_ref_reps_arraysize_q = q{CREATE INDEX IF NOT EXISTS
        "idx_fasta_ref_reps_arraysize" ON
        "fasta_ref_reps" (lastindex - firstindex + 1)};

    # $dbh->do($create_fasta_ref_reps_q);
    # $dbh->do($fasta_ref_reps_index_q);
    # $dbh->do($fasta_ref_reps_patsize_q);
    # $dbh->do($fasta_ref_reps_arraysize_q);

    # If the table does not exist, or we are forcing a redo, load
    # the profiles into the db.
    my $tab_exists_sth
        = $dbh->table_info( undef, "main", "fasta^_ref^_reps", "TABLE",
        { Escape => '^' } );

    my $tab_exists
        = scalar @{ $tab_exists_sth->fetchall_arrayref( undef, 1 ) };

    ( $ENV{DEBUG} )
        && warn "Does ref table exist?: $tab_exists. Redo set to: $redo\n";

    if ( $tab_exists == 0 || $redo ) {
        print "Creating reference sequence database.\n";
        $dbh->do(q{DROP TABLE IF EXISTS fasta_ref_reps});
        $dbh->do($create_fasta_ref_reps_q);
        $dbh->do($fasta_ref_reps_index_q);
        $dbh->do($fasta_ref_reps_patsize_q);
        $dbh->do($fasta_ref_reps_arraysize_q);

        # Read in ref file and create a virtual table out of it
        $dbh->sqlite_create_module(
            perl => "DBD::SQLite::VirtualTable::PerlData" );

        our $seq_rows = [];

        # my $installdir = "$FindBin::RealBin";
        open my $refset, "<", $refseq
            or confess "Error opening reference sequences file "
            . $refseq
            . ": $!.\nStopped at";

        # Skip header
        <$refset>;

        # Insert all ref TRs from refset file
        while ( my $r = <$refset> ) {
            chomp $r;
            my @fields = split /,/, $r;
            push @$seq_rows, \@fields;
        }

        ( $ENV{DEBUG} )
            && warn "Read "
            . scalar(@$seq_rows)
            . " lines from refseq file\n";

        close $refset;

        # Create a virtual tables for the inputs
        $dbh->begin_work;
        $dbh->do(
            q{CREATE VIRTUAL TABLE temp.seqtab
            USING perl(rid integer,
                firstindex integer,
                lastindex integer,
                copynum float,
                head varchar(100),
                flankleft text,
                pattern text,
                sequence text,
                flankright text,
                conserved float,
            arrayrefs="vutil::seq_rows")}
        );
        $dbh->commit;

        my $cols = join(
            ",", qw(rid
                firstindex
                lastindex
                copynum
                head
                flankleft
                pattern
                sequence
                flankright
                conserved)
        );

        $dbh->begin_work;
        my $num_rows = $dbh->do(
            qq{INSERT INTO fasta_ref_reps ($cols)
            SELECT * FROM temp.seqtab}
        );

        if ( $num_rows != @$seq_rows ) {
            $dbh->rollback;

            # unlink($refdbfile);
            die "Error inserting reference sequences into "
                . $dbh->sqlite_db_filename()
                . ": inserted $num_rows but read "
                . scalar(@$seq_rows)
                . " lines from file. Aborting...";

            $dbh->disconnect;
        }

        $dbh->commit;
    }
}

################################################################
# Use this to load reference set db with profiles (leb36 file)
sub load_refprofiles_db {
    my ( $dbh, $refleb36, $redo ) = @_;

    # By default, set redund flag to 1 so we can set
    # non-redundant trs to 0 later;
    my $create_ref_profiles_q = q{CREATE TABLE IF NOT EXISTS `ref_profiles` (
    `rid` integer NOT NULL,
    `proflen` integer NOT NULL,
    `proflenrc` integer NOT NULL,
    `profile` text NOT NULL,
    `profilerc` text NOT NULL,
    `nA` integer NOT NULL,
    `nC` integer NOT NULL,
    `nG` integer NOT NULL,
    `nT` integer NOT NULL,
    `redund` integer NOT NULL default 1,
    `minrepidx` integer default 0,
    PRIMARY KEY (`rid`))};
    $dbh->do($create_ref_profiles_q);

    # If the table does not exist, or we are forcing a redo, load
    # the profiles into the db.
    my ($num_rows)
        = $dbh->selectrow_array(q{SELECT COUNT(*) FROM ref_profiles});
    if ( $num_rows == 0 || $redo ) {
        $dbh->begin_work;
        $dbh->do(q{DROP TABLE IF EXISTS `ref_profiles`});
        $dbh->do($create_ref_profiles_q);
        $dbh->commit;

        # Read in ref leb36 file and create a virtual table out of it
        try {
            $dbh->sqlite_create_module(
                perl => "DBD::SQLite::VirtualTable::PerlData" );
        }
        catch {
            if (/sqlite_create_module failed with error not an error/) {
                ( $ENV{DEBUG} )
                    && warn
                    "Not creating VirtualTable; module already registered.\n";
            }
            else {
                die
                    "Error installing VirtualTable in SQLite db handle: $_\n.";
            }
        };

        our $leb36_rows = [];
        open my $refleb36_fh, "<", "$refleb36"
            or die "Error opening reference profiles file: $!\n";

        while ( my $r = <$refleb36_fh> ) {
            chomp $r;
            my @fields = split /\s+/, $r;
            push @$leb36_rows, [ @fields[ 0, 3 .. 11 ] ];
        }

        ( $ENV{DEBUG} )
            && warn "Read "
            . scalar(@$leb36_rows)
            . " lines from refleb file\n";

        close $refleb36_fh;

        $dbh->begin_work;
        $dbh->do(
            q{CREATE VIRTUAL TABLE temp.leb36tab
            USING perl(rid integer,
                proflen integer,
                proflenrc integer,
                profile text,
                profilerc text,
                nA integer,
                nC integer,
                nG integer,
                nT integer,
            arrayrefs="vutil::leb36_rows")}
        );
        $dbh->commit;

        my $cols = join(
            ",", qw(rid
                proflen
                proflenrc
                profile
                profilerc
                nA
                nC
                nG
                nT)
        );

        # Insert all ref TRs from refleb36 file
        $dbh->begin_work();
        my $num_rows = $dbh->do(
            qq{INSERT INTO ref_profiles ($cols)
            SELECT * FROM temp.leb36tab}
        );

        if ( $num_rows != @$leb36_rows ) {
            $dbh->rollback();
            $dbh->disconnect();
            die
                "Error inserting reference profiles into database: inserted $num_rows but read "
                . scalar(@$leb36_rows)
                . " lines from file. Aborting...";
        }

        $dbh->do(q{DROP TABLE temp.leb36tab});

        $dbh->commit();

        # Run redund
        try {
            print "Running redund on ref set.\n";
            run_redund( $dbh, $refleb36, "reference.leb36", 0, 0 );
        }
        catch {
            warn "Error running redund on reference set: $_\n";
            # if (/DBD::SQLite::db/) {
            # }

            $dbh->rollback();
            $dbh->do(q{DROP TABLE IF EXISTS `ref_profiles`});
            $dbh->disconnect();

            die "\n";
        };
        return 1;
    }

    return 0;
}

################################################################
sub set_indist {
    my ( $dbh, $indist_file, $redo ) = @_;
    my ($indists_unset) = $dbh->selectrow_array(
        q{SELECT COUNT(*) IN
    (SELECT COUNT(*) FROM fasta_ref_reps
    GROUP BY is_indist)
    FROM fasta_ref_reps}
    );

    return
        unless ( $redo || $indists_unset );

    print "Updating indistinguishable TRs.\n";
    open my $indist_fh, "<", $indist_file
        or die "Error opening indist file $indist_file: $!\n";
    my $count = 0;
    while ( my $line = <$indist_fh> ) {
        $count++;
        $line =~ /-(\d+)/;

        $dbh->do(
            qq{UPDATE fasta_ref_reps SET is_singleton=0,is_indist=1 WHERE rid=$1}
        ) or die $dbh->errstr;
    }
}

################################################################
# Use database to produce the profiles needed to run redund.exe
# and update database on redundant TRs.
sub run_redund {

    # First, create a temporary directory and copy the filtered set there
    my $tmpdir      = File::Temp->newdir();
    my $tmpdir_name = $tmpdir->dirname;
    my ( $dbh, $input_refset, $output_basename, $keep_files, $redo ) = @_;

    # # Add - sign if negating for ref set and set new input to same
    # # file path as the final output
    # if ($use_as_ref) {
    #     my $cmd = q(awk '{print "-"$0}')
    #         . qq{ $input_refset > $tmpdir_name/$output_basename};

    #     # warn "awk command: $cmd\n";
    #     system($cmd);
    #     $input_refset = "$tmpdir_name/$output_basename";

    #     # warn "New input: $input_refset\n";
    # }

    # Check if we need to run redund, or are forcing a rerun

    my ($num_rows)
        = $dbh->selectrow_array(
        q{SELECT COUNT(*) FROM ref_profiles WHERE redund = 0});

    # Return if redund has already been run for this refset
    # and were are not forcing a redo
    return
        unless ( $num_rows == 0 || $redo );

    #=<<Run redund.exe on the input leb36 files.>>
    # First run sorts by minimum representation
    my $tmp_file = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );
    my $installdir = "$FindBin::RealBin";
    system("$installdir/redund.exe $input_refset $tmp_file -s -i");
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    # Second run eliminates redundancies
    system(
        "$installdir/redund.exe $tmp_file $tmpdir_name/$output_basename -n -i"
    );
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    #=<<Mark all TRs in file as non-redundant>>
    my $rotindex_fn = "$tmpdir_name/$output_basename.rotindex";
    open my $set_fh, "<", $rotindex_fn;

    my @unique_trs;
    while ( my $entry = <$set_fh> ) {
        $entry =~ /(\d+)['"]/;
        my $trid = $1;
        push @unique_trs, $trid;
    }
    close $set_fh;

    # warn Dumper($unique_trs) ."\n";

    print "Marking non-redundant TRs in database.\n";
    $dbh->begin_work();
    $dbh->do(
        q{CREATE TEMPORARY TABLE temp.unique_trs
        (`rid` integer PRIMARY KEY)}
    ) or die "Couldn't do statement: " . $dbh->errstr;

    my $insert_sth
        = $dbh->prepare(q{INSERT INTO temp.unique_trs (rid) VALUES (?)});
    for my $rid (@unique_trs) {
        $insert_sth->execute($rid);
    }

    # Update the ref_profiles table
    my $update_redund = q{UPDATE ref_profiles SET redund=0
        WHERE EXISTS (
        SELECT rid FROM temp.unique_trs t2
        WHERE ref_profiles.rid = t2.rid
    )};
    my $unique_tr_count = $dbh->do($update_redund)
        or die "Couldn't do statement: " . $dbh->errstr;
    $dbh->commit;

    $dbh->do(
        q{CREATE TABLE IF NOT EXISTS `files` (
        `rotindex` BLOB
        )}
    );

    #=<<Set minimum representation>>
    # Attach database created by redund.exe
    print "Getting sort order of TRs by minimum representation.\n";
    my $minreporder_db = "$tmp_file.db";
    $dbh->do(qq{ATTACH DATABASE "$minreporder_db" AS minrep})
        or carp "Couldn't do statement: $DBI::errstr\n";
    ( $ENV{DEBUG} ) && warn "Connecting to SQLite db at $minreporder_db\n";

    # print "Press enter to continue...";
    # my $dummy = <STDIN>;

    $dbh->begin_work();

    # Copy the CREATE TABLE query for the minreporder table
    # and create the same table in the ref set db, copying
    # over the data
    my ($minrep_sql) = $dbh->selectrow_array(
        q{SELECT sql FROM minrep.sqlite_master
        WHERE name = 'minreporder'}
    );

    unless ($minrep_sql) {
        my $res = $dbh->selectall_arrayref(q{SELECT * FROM minrep.sqlite_master});
        use Data::Dumper;
        print Dumper($res) . "\n";
        die "Unable to get CREATE TABLE statement for minreporder from redund db\n";
    }

    $minrep_sql =~ s/minreporder/main.minreporder/;
    $dbh->do(q{DROP TABLE IF EXISTS main.minreporder})
        or carp "Couldn't do statement: $DBI::errstr\n";
    $dbh->do($minrep_sql)
        or carp "Couldn't do statement: $DBI::errstr\n";
    $dbh->do(
        q{INSERT INTO main.minreporder
            SELECT * FROM minrep.minreporder}
    ) or carp "Couldn't do statement: $DBI::errstr\n";
    $dbh->commit;
    $dbh->do(qq{DETACH minrep})
        or carp "Couldn't do statement: $DBI::errstr\n";

    # Save the rotindex file
    print "Saving rotindex (redundancy index) into database.\n";
    my $rotindex;
    {
        local $/;
        open my $fh, '<', $rotindex_fn or die "can't open $rotindex_fn: $!";
        $rotindex = <$fh>;
        close $fh;
    }
    my $load_rotindex_sth
        = $dbh->prepare(qq{INSERT INTO files (rotindex) VALUES (?)})
        or die "Couldn't do statement: $DBI::errstr";
    $load_rotindex_sth->bind_param( 1, $rotindex, DBI::SQL_BLOB );
    $load_rotindex_sth->execute;

    print "$unique_tr_count unique TRs in filtered set.\n";

    unlink($minreporder_db);

    if ($keep_files) {
        return $tmpdir;
    }
}

################################################################
# TODO Function to generate leb36 file from db after redund was
# run. Will execute run_redund unless 'minreporder' table exists
# SQL query: SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name='minreporder';
# Code to import from produce_indist

################################################################
sub write_sqlite {
    my $installdir = "$FindBin::RealBin";
    ( $ENV{DEBUG} ) && warn "Creating SQLite database...\n";

    my $dbfile = "$cnf{OUTPUT_DIR}/vntr_$cnf{RUN_NAME}/$cnf{RUN_NAME}.db";
    system("sqlite3 $dbfile < $installdir/sqlite_schema.sql");
}

################################################################
sub sqlite_install_RC_function {
    my ($dbh) = @_;
    $dbh->sqlite_create_function(
        'RC', 1,
        sub {
            # complement reversed DNA sequence
            my $seq = shift;

            $seq = reverse $seq;

            $seq =~ tr/ACGT/TGCA/;

            return $seq;
        }
    );
}

################################################################
sub gen_exec_array_cb {
    my $arr_ref = shift
        or croak "Error: Must specify an array reference.\n";

    my $arr_idx = 0;
    return sub {
        return ( $arr_idx < @$arr_ref ) ? $arr_ref->[ $arr_idx++ ] : undef;
    }
}

#-------------------------------------------------------------------------------
## @brief      vs_db_insert
##
## @param      $dbh A database handle
## @param      $sth A statement handle for the INSERT statement you
##             wish to perform.
## @param      $cb_or_sth Either: a callback function which returns
##             an arrayref representing a row to insert on each call
##             OR, a statement handle which will return rows to be
##             inserted by $sth into the database.
## @param      $errstr A message to print in case of error.
##
## @return     Upon success, the number of rows inserted into the
##             database.
##
sub vs_db_insert {
    my ( $dbh, $sth, $cb_or_sth, $errstr ) = @_;
    my ( $tuples, @tuple_status );

    $dbh->begin_work();
    try {
        $tuples = $sth->execute_array(
            {   ArrayTupleFetch  => $cb_or_sth,
                ArrayTupleStatus => \@tuple_status
            }
        );
        ( $ENV{DEBUG} ) && warn "Inserted $tuples records into database\n";
    }
    catch {
        warn "$_\n";

        use Set::IntSpan;
        my %status_hash;

        for my $tuple ( 0 .. @tuple_status - 1 ) {
            my $status = $tuple_status[$tuple];
            $status = [ 0, "Skipped" ]
                unless defined $status;
            next unless ref $status;
            push @{ $status_hash{ $status->[1] } }, $tuple;
        }

        while ( my ( $msg, $rows ) = each %status_hash ) {
            warn "$msg (rows: "
                . Set::IntSpan->new($rows)->run_list() . ")\n";
        }

        eval { $dbh->rollback(); };
        if ($@) {
            warn "Database rollback failed.\n";
        }
    };

    croak "$errstr\n"
        unless ($tuples);
    $dbh->commit();
    return $tuples;
}

################################################################



####################################################################################

1;

