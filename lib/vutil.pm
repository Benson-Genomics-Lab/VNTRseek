package vutil;
use v5.24;
use Cwd;
use DBI;
use Carp;
use FindBin;
use File::Temp;
use Try::Tiny;
use POSIX qw(strftime);

if ( $ENV{DEBUG} ) {
    use Data::Dumper;
}

use base 'Exporter';
our @EXPORT_OK
    = qw(read_config_file get_config set_config set_credentials get_dbh get_ref_dbh make_refseq_db load_refprofiles_db run_redund write_sqlite set_statistics get_statistics set_datetime print_config trim create_blank_file get_trunc_query sqlite_install_RC_function gen_exec_array_cb vs_db_insert);

# vutil.pm
# author: Yevgeniy Gelfand, Yozen Hernandez
# create date: Oct 30, 2010
# function: common functions for vntr pipleline scripts.
# Provide functions for database management, configuration
# reading/writing, and miscellaneous common utilities.

my %VSCNF_FILE  = ();
my $VSREAD      = 0;
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
    ERROR_STEP                                          => 1,
    ERROR_DESC                                          => 1,
    ERROR_CODE                                          => 1,
    N_MIN_SUPPORT                                       => 1,
    MIN_FLANK_REQUIRED                                  => 1,
    MAX_FLANK_CONSIDERED                                => 1,
);

################################################################
sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub create_blank_file {
    my $filename = shift or croak("Error: filename is a required argument\n");
    open my $blank_file, ">", $filename
        or croak("Error creating blank file $filename: $!");
    close $blank_file or croak("Error closing blank file $filename: $!");
}

################################################################
sub read_config_file {

    # Get file location
    unless (@_) { die "read_config_file: expects 1 parameters.\n"; }
    my $file_loc = shift;

    ( $ENV{DEBUG} ) && warn "Reading configuration file: $file_loc\n";

    # read global config
    my $cnf;
    try {
        open( $cnf, "<", "$file_loc" )
            or die "$!";
    }
    catch {
        ( $ENV{DEBUG} ) && warn "$_\n";
        $VSCNF_FILE{ERR} = $_;
        return 0;
    };

    while (<$cnf>) {
        chomp;

        # skip start comments
        if (/^\#/) { next; }

        if (/(.+)=(.*)/) {
            my $key = trim($1);
            my $val = trim($2);
            $val =~ s/\s*\#.*//;    # strip end comments
            $VSCNF_FILE{ uc($key) } = $val;

            #print $key."=".$val."\n";
        }
    }

    close($cnf);
    return 1;

}

################################################################
sub get_config {
    croak "Error: function expects 2 parameters\n"
        unless ( @_ == 2 );

    my ( $dbsuffix, $config_loc ) = @_;
    my $installdir  = "$FindBin::RealBin";
    my $config_file = "$config_loc/$dbsuffix.vs.cnf";
    $VSCNF_FILE{ERR} = "";

    # $VSCNF_FILE{NEW_RUN} = 0;
    unless ($VSREAD) {
        $VSCNF_FILE{ERR} = "";

        # Must read global file first. Sets up the defaults.
        read_config_file("$installdir/vs.cnf");
        if ( $VSCNF_FILE{ERR} =~ /Permission denied/ ) {
            die "Could not read global config. Please make sure "
                . "your user or group has read permissions.\n";
        }
        elsif ( $VSCNF_FILE{ERR} =~ /No such file or directory/ ) {
            die "Global config does not exist! Either reinstall "
                . "VNTRseek or copy the default configuration from "
                . "source distribution and modify as needed.\n";
        }
        elsif ( $VSCNF_FILE{ERR} ) {
            die "Unexpected error reading global configuration file: "
                . $VSCNF_FILE{ERR} . "\n";
        }

        $VSCNF_FILE{ERR} = "";
        read_config_file($config_file);
        if ( $VSCNF_FILE{ERR} =~ /No such file or directory/ ) {
            warn "Run config does not exist. "
                . "A new one will be created, but make sure your "
                . "run name (dbsuffix) is correct!\n";

            # $VSCNF_FILE{NEW_RUN} = 1;
        }
        elsif ( $VSCNF_FILE{ERR} =~ /Permission denied/ ) {
            die "Could not read run config. Please make sure "
                . "your user or group has read AND write permissions.\n";
        }
        elsif ( $VSCNF_FILE{ERR} ) {
            die "Unexpected error reading run configuration file: "
                . $VSCNF_FILE{ERR} . "\n";
        }
        delete $VSCNF_FILE{ERR};

        # Set VSREAD;
        $VSREAD               = 1;
        $VSCNF_FILE{DBSUFFIX} = $dbsuffix;
        $VSCNF_FILE{CONF_DIR} = $config_loc;
    }

    return %VSCNF_FILE;
}

################################################################
sub set_config {
    my %in_hash = @_;

    if ( $ENV{DEBUG} ) {
        use Data::Dumper;
        warn Dumper( \%in_hash );
    }

    # Validation
    unless ( $in_hash{SERVER} ) {
        croak(
            "Please set machine name (SERVER) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{NPROCESSES} <= 0 ) {
        warn(     "NPROCESSES not set. Using default of 2. "
                . "Please set number of processes to be used by the pipeline "
                . "(NPROCESSES) variable on the command line or in the "
                . "configuration file.\n " );
        $in_hash{NPROCESSES} = 2;
    }

    if ( $in_hash{MIN_FLANK_REQUIRED} <= 0 ) {
        croak(
            "Please set min flank required to be used by the pipeline (MIN_FLANK_REQUIRED) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{MAX_FLANK_CONSIDERED} <= 0 ) {
        croak(
            "Please set max flank required to be used by the pipeline (MAX_FLANK_CONSIDERED) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{MIN_SUPPORT_REQUIRED} <= 0 ) {
        croak(
            "Please set min support required to be used by the pipeline (MIN_SUPPORT_REQUIRED) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{STRIP_454_KEYTAGS} < 0 ) {
        croak(
            "Please set strip_454_keytags (strip_454_keytags) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{REFERENCE_INDIST_PRODUCE} < 0 ) {
        croak(
            "Please set reference_indist_produce to be used by the pipeline (reference_indist_produce) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{IS_PAIRED_READS} < 0 ) {
        croak(
            "Please set is_paired_reads (is_paired_reads) variable on the command line or in the configuration file.\n"
        );
    }

    unless ( $in_hash{PLOIDY} ) {
        $in_hash{PLOIDY} = 2;
        warn(
            "Option 'ploidy' is not set. Setting default: $in_hash{PLOIDY}.\n"
        );
    }

    if ( !defined $in_hash{HTML_DIR} || ( "" eq $in_hash{HTML_DIR} ) ) {

        # warn "Warning: 'html_dir' option is unset.\n";
    }
    elsif ( !( -e $in_hash{HTML_DIR} ) && !mkdir("$in_hash{HTML_DIR}") ) {
        croak("Could not create html_dir directory ($in_hash{HTML_DIR}).\n");
    }

    # For older runs where FASTA_DIR was the name of the option
    if ( exists $in_hash{FASTA_DIR} && $in_hash{FASTA_DIR} ) {
        $in_hash{INPUT_DIR} = $in_hash{FASTA_DIR};
    }

    unless ( exists $in_hash{INPUT_DIR} && $in_hash{INPUT_DIR} ) {
        croak(
            "Please set input directory (input_dir) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{OUTPUT_ROOT} ) {
        croak(
            "Please set output directory (output_root) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{REFERENCE} ) {
        croak(
            "Please set reference file base name (reference) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{TMPDIR} ) {
        warn
            "'tmpdir' is not set. Using '/tmp'. Call vntrseek with --tmpdir <directory> to change this.\n";
    }

    if ( !( -e $in_hash{REFERENCE} . ".db" ) || $in_hash{REDO_REFDB} ) {
        my $err;
        $in_hash{REFERENCE_DB}     = $in_hash{REFERENCE} . ".db";
        $in_hash{REFERENCE_SEQ}    = $in_hash{REFERENCE} . ".seq";
        $in_hash{REFERENCE_FILE}   = $in_hash{REFERENCE} . ".leb36";
        $in_hash{REFERENCE_INDIST} = $in_hash{REFERENCE} . ".indist";
        unless ( -e $in_hash{REFERENCE_FILE} ) {
            $err = "File '$in_hash{REFERENCE_FILE}' not found!";
        }
        unless ( -e $in_hash{REFERENCE_SEQ} ) {
            $err = "File '$in_hash{REFERENCE_SEQ}' not found!";
        }
        if (   !( -e $in_hash{REFERENCE_INDIST} )
            && !$in_hash{REFERENCE_INDIST_PRODUCE} )
        {
            $err = "File '$in_hash{REFERENCE_INDIST}' not found!";
        }

        die "Error: $err File necessary to build reference set database "
            . "(see VNTRseek documentation).\n"
            if ($err);
    }

    unless ( -e $in_hash{INPUT_DIR} ) {
        die("Error: input directory '$in_hash{INPUT_DIR}' not found!");
    }
    if ( !-e $in_hash{OUTPUT_ROOT} && !mkdir( $in_hash{OUTPUT_ROOT} ) ) {
        die("Error creating output root ($in_hash{OUTPUT_ROOT}). Please check the path or manually create this directory and try again."
        );
    }
    unless ( -e $in_hash{OUTPUT_ROOT} ) {
        die("Directory '$in_hash{OUTPUT_ROOT}' not found!");
    }
    unless ( -e $in_hash{TMPDIR} ) {
        die("Temporary directory '$in_hash{TMPDIR}' not found!");
    }
    if ( $in_hash{HTML_DIR} && !( -x $in_hash{HTML_DIR} ) ) {
        die("Directory '$in_hash{HTML_DIR}' not executable!");
    }
    unless ( -x $in_hash{OUTPUT_ROOT} ) {
        die("Directory '$in_hash{OUTPUT_ROOT}' not executable!");
    }
    unless ( -x $in_hash{TMPDIR} ) {
        die("Directory '$in_hash{TMPDIR}' not executable!");
    }

    if ( $in_hash{HTML_DIR} && !( -w $in_hash{HTML_DIR} ) ) {
        die("Directory '$in_hash{HTML_DIR}' not writable!");
    }
    unless ( -w $in_hash{OUTPUT_ROOT} ) {
        die("Directory '$in_hash{OUTPUT_ROOT}' not writable!");
    }
    unless ( -w $in_hash{TMPDIR} ) {
        die("Directory '$in_hash{TMPDIR}' not writable!");
    }

    # DELETEME While testing, this option is always sqlite.
    # In final version, this option will be gone completely.
    $in_hash{BACKEND} = "sqlite";
    %VSCNF_FILE       = %in_hash;
    $VSREAD           = 1;
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
    unless ($VSREAD) {
        carp
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    # my $DBSUFFIX = shift;
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

    if ( $opts->{redo} && $VSREAD ) {
        $VSCNF_FILE{REDO_REFDB} = 0;
        print_config( $VSCNF_FILE{CONF_DIR} );
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
    unless ( @_ >= 1 ) {
        croak "Error: expecting at least 1 argument, got none.\n";
    }
    my ( $refbase, $opts ) = @_;
    ( $ENV{DEBUG} ) && warn "get_ref_dbh arguments: " . Dumper( \@_ ) . "\n";
    my $dbh = _init_ref_dbh( $refbase, $opts );

    return $dbh;
}

################################################################
sub get_dbh {
    unless ($VSREAD) {
        croak
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    # Might use in the future, connection options
    my $opts = shift // {};

    my $dbh;
    my $dbfile
        = "$VSCNF_FILE{OUTPUT_ROOT}/vntr_$VSCNF_FILE{DBSUFFIX}/$VSCNF_FILE{DBSUFFIX}.db";

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

    $dbh
        = DBI->connect( "DBI:SQLite:dbname=$dbfile",
        undef, undef, \%dbi_opts, )
        or die "Could not connect to database $dbfile: $DBI::errstr";

    my ($stats_schema) = $dbh->selectrow_array(
        q{SELECT sql FROM sqlite_master
        WHERE name = 'stats'}
    );

    # If the stats table does not exist, init the db.
    unless ($stats_schema) {
        write_sqlite();
    }

    my ($schema_ver) = $dbh->selectrow_array(q{PRAGMA user_version});

    if ( $schema_ver < 1 ) {
        $dbh->disconnect;
        $dbh
            = DBI->connect( "DBI:SQLite:dbname=$dbfile", undef, undef,
            { %dbi_opts, ReadOnly => 0, } )
            or die "Could not connect to database $dbfile: $DBI::errstr";
        $dbh->do(q{PRAGMA foreign_keys = off});
        $dbh->begin_work;
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
        $dbh->commit;
        $dbh->do(q{PRAGMA foreign_keys = on});
        $dbh->do(q{PRAGMA user_version = 1});
        $schema_ver = 1;
        $dbh->disconnect;
        $dbh
            = DBI->connect( "DBI:SQLite:dbname=$dbfile",
            undef, undef, \%dbi_opts, )
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
        my $refdbfile = $VSCNF_FILE{REFERENCE} . ".db";

        # First initialize the reference DB, if needed
        # Don't save the returned dbh
        # _init_ref_dbh( $VSCNF_FILE{REFERENCE},
        #     { redo => $VSCNF_FILE{REDO_REFDB} } );
        ( $ENV{DEBUG} )
            && warn "Connecting to refseq db at $refdbfile\n";
        try {
            $dbh->do(qq{ATTACH DATABASE "$refdbfile" AS refdb});
        }
        catch {
            warn "$_\n";
            if (/unable to open database/) {
                die "Try running VNTRseek as:\n"
                    . "\n\tvntrseek --reference $VSCNF_FILE{REFERENCE} --redo_refdb\n\n"
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
        warn "Creating reference sequence database...\n";
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

####################################
sub get_trunc_query {
    die "Error: need table name for truncate query\n"
        unless @_ == 2;
    my ( $backend, $table ) = @_;
    my $trunc_query;
    if ( $backend eq "sqlite" ) {
        $trunc_query = qq{DELETE FROM $table};
    }
    elsif ( $backend eq "mysql" ) {
        $trunc_query = qq{TRUNCATE TABLE $table};
    }

    return $trunc_query;
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
        $dbh->begin_work;
        my $num_rows = $dbh->do(
            qq{INSERT INTO ref_profiles ($cols)
            SELECT * FROM temp.leb36tab}
        );

        if ( $num_rows != @$leb36_rows ) {
            $dbh->rollback;
            $dbh->disconnect;
            die
                "Error inserting reference profiles into database: inserted $num_rows but read "
                . scalar(@$leb36_rows)
                . " lines from file. Aborting...";
        }

        $dbh->do(q{DROP TABLE temp.leb36tab});

        $dbh->commit;

        # Run redund
        try {
            warn "Running redund on ref set\n";
            run_redund( $dbh, $refleb36, "reference.leb36", 0, 0 );
        }
        catch {
            warn "Error running redund on reference set: $_\n";
            # if (/DBD::SQLite::db/) {
            # }

            $dbh->rollback;
            $dbh->do(q{DROP TABLE IF EXISTS `ref_profiles`});
            $dbh->disconnect;

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

    warn "Updating indistinguishable TRs...\n";
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

    warn "Marking non-redundant TRs in database.\n";
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
    warn "Getting sort order of TRs by minimum representation.\n";
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
        say Dumper($res);
        die "unable to get CREATE TABLE statement for minreporder from redund db\n";
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
    warn "Saving rotindex (redundancy index) into database\n";
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

    warn "$unique_tr_count unique TRs in filtered set\n";

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
    unless ($VSREAD) {
        carp
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    my $output_folder = "$VSCNF_FILE{OUTPUT_ROOT}/vntr_$VSCNF_FILE{DBSUFFIX}";
    if ( !-e "$output_folder" ) {
        unless ( mkdir("$output_folder") ) {
            if ( $!{EEXIST} ) {
                warn
                    "Warning: Not creating data directory, directory exists.\n";
            }
            else {
                warn "Warning: Failed to create data directory!\n";
            }
        }
    }

    my $installdir = "$FindBin::RealBin";
    ( $ENV{DEBUG} ) && warn "Creating SQLite database...\n";

    my $dbfile
        = "$VSCNF_FILE{OUTPUT_ROOT}/vntr_$VSCNF_FILE{DBSUFFIX}/$VSCNF_FILE{DBSUFFIX}.db";

# if ( -e $dbfile && !( exists $VSCNF_FILE{CLEAN} ) ) {
#     die
#         "Error: run database exists. Running step 0 is potentially destructive! "
#         . "Exiting...\n(To override use option '--clean')\n";
# }
# ( exists $VSCNF_FILE{CLEAN} ) && unlink $dbfile;
    my $exestring = "sqlite3 $dbfile < $installdir/sqlite_schema.sql";
    warn "Executing: $exestring\n";
    system($exestring);
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
        or croak "Error: Must specify an array reference. "
        . "(Developer made a mistake here)";

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

    $dbh->begin_work;
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

        eval { $dbh->rollback; };
        if ($@) {
            warn "Database rollback failed.\n";
        }
    };

    croak "$errstr\n"
        unless ($tuples);
    $dbh->commit;
    return $tuples;
}

################################################################

sub print_config {

    my $argc = @_;
    if ( $argc < 1 ) {
        croak "print_config: expects 1 parameters, passed $argc!\n";
    }

    my $startdir = $_[0];

    if ( !defined $VSCNF_FILE{"NPROCESSES"} ) {
        $VSCNF_FILE{"NPROCESSES"} = -1;
    }
    if ( !defined $VSCNF_FILE{"MIN_FLANK_REQUIRED"} ) {
        $VSCNF_FILE{"MIN_FLANK_REQUIRED"} = -1;
    }
    if ( !defined $VSCNF_FILE{"MAX_FLANK_CONSIDERED"} ) {
        $VSCNF_FILE{"MAX_FLANK_CONSIDERED"} = -1;
    }
    if ( !defined $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"} ) {
        $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"} = -1;
    }

    if ( !defined $VSCNF_FILE{"SERVER"} ) { $VSCNF_FILE{"SERVER"} = ""; }
    if ( !defined $VSCNF_FILE{"STRIP_454_KEYTAGS"} ) {
        $VSCNF_FILE{"STRIP_454_KEYTAGS"} = -1;
    }
    if ( !defined $VSCNF_FILE{"IS_PAIRED_READS"} ) {
        $VSCNF_FILE{"IS_PAIRED_READS"} = -1;
    }
    if ( !defined $VSCNF_FILE{"HTML_DIR"} ) { $VSCNF_FILE{"HTML_DIR"} = ""; }
    if ( !defined $VSCNF_FILE{"INPUT_DIR"} ) {
        $VSCNF_FILE{"INPUT_DIR"} = "";
    }
    if ( !defined $VSCNF_FILE{"OUTPUT_ROOT"} ) {
        $VSCNF_FILE{"OUTPUT_ROOT"} = "";
    }
    if ( !defined $VSCNF_FILE{"TMPDIR"} ) { $VSCNF_FILE{"TMPDIR"} = ""; }
    if ( !defined $VSCNF_FILE{"REFERENCE_FILE"} ) {
        $VSCNF_FILE{"REFERENCE_FILE"} = "";
    }
    if ( !defined $VSCNF_FILE{"REFERENCE_SEQ"} ) {
        $VSCNF_FILE{"REFERENCE_SEQ"} = "";
    }
    if ( !defined $VSCNF_FILE{"REFERENCE_INDIST"} ) {
        $VSCNF_FILE{"REFERENCE_INDIST"} = "";
    }
    if ( !defined $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"} ) {
        $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"} = -1;
    }

    # look in the directory the script was started in
    my $config_file = "$startdir/" . $VSCNF_FILE{DBSUFFIX} . ".vs.cnf";
    if ( open( my $config_fh, ">", $config_file ) ) {

        print $config_fh <<CNF;
# Database backend
BACKEND=$VSCNF_FILE{"BACKEND"}

# set this to the number of processors on your system 
# (or less if sharing the system with others or RAM is limited)
# eg, 8
NPROCESSES=$VSCNF_FILE{"NPROCESSES"}

# minimum required flank on both sides for a read TR to be considered
# eg, 10
MIN_FLANK_REQUIRED=$VSCNF_FILE{"MIN_FLANK_REQUIRED"}

# maximum flank length used in flank alignments
# set to big number to use all
# if read flanks are long with a lot of errors, 
# it might be useful to set this to something like 50
# max number of errors per flank is currently set to 8 (can be changed in main script only)
# eg, 1000
MAX_FLANK_CONSIDERED=$VSCNF_FILE{"MAX_FLANK_CONSIDERED"}

# minimum number of mapped reads which agree on copy number to call an allele
# eg, 2
MIN_SUPPORT_REQUIRED=$VSCNF_FILE{"MIN_SUPPORT_REQUIRED"}

# Whether or not to keep reads detected as PCR duplicates. A nonzero (true) value
# means that detected PCR duplicates will not be removed. Default is 0.
KEEPPCRDUPS=$VSCNF_FILE{"KEEPPCRDUPS"}

# server name, used for html generating links
# eg, orca.bu.edu
SERVER=$VSCNF_FILE{"SERVER"}

# for 454 platform, strip leading 'TCAG' 
# eg, 1 - yes
# eg, 0 - no (use no for all other platforms)
STRIP_454_KEYTAGS=$VSCNF_FILE{"STRIP_454_KEYTAGS"}

# data is paired reads
# eg, 0 = no 
# eg, 1 - yes
IS_PAIRED_READS=$VSCNF_FILE{"IS_PAIRED_READS"}

# Sample ploidy. Default is 2. For haploid, set to 1.
PLOIDY=$VSCNF_FILE{"PLOIDY"}

# Rebuild reference database
# eg, 0 = no 
# eg, 1 - yes
REDO_REFDB=$VSCNF_FILE{"REDO_REFDB"}

# input data directory 
# (plain or gzipped fasta/fastq files)
# eg, /input
INPUT_DIR=$VSCNF_FILE{"INPUT_DIR"}

# output directory (must be writable and executable!)
# eg, /output
OUTPUT_ROOT=$VSCNF_FILE{"OUTPUT_ROOT"}

# temp (scratch) directory (must be executable!)
# eg, /tmp
TMPDIR=$VSCNF_FILE{"TMPDIR"}

# names for the reference files 

# (leb36 file, sequence plus flank data file, indistinguishable references file) 
# files must be in install directory

# eg, hg19. This is the base name for files describing
# reference TR loci (.db, .seq, .leb36, and .indist)
REFERENCE=$VSCNF_FILE{"REFERENCE"}

# generate a file of indistinguishable references, 
# necessary only if a file is not already available for the reference set
# eg, 1- generate
# eg, 0 - don't generate
REFERENCE_INDIST_PRODUCE=$VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"}

CNF

        close($config_fh);

        chmod 0640, $config_file;

    }
    else {

        die "print_config: can't open '$config_file' for writing!\n";
    }

}

####################################################################################

1;

