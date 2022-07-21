#!/usr/bin/env perl

# command line usage example:
#  ./run_proclu.pl 1 inputfolder referencefolder " 88 " 8 psearch.exe 4
# where 1 is the number of files to process in one batch
# and inputfolder is the input directory containing leb36 files
# note: cluster params is a hard cutoff QUOTES ARE REQUIRED
# 8 here is number of processes
# psearch.exe is name of clustering executable
# 4 is maxerror in either flank

use strict;
use warnings;
use Cwd;
use POSIX qw(strftime);

warn strftime( "Start: %F %T\n\n", localtime );

die "Useage: run_proclu.pl expects 8 arguments.\n"
    unless scalar @ARGV >= 8;

my $curdir           = getcwd();
my $files_to_process = $ARGV[0]; # number of files to process in one batch
my $tgz_dir          = $ARGV[1];
my $reffolder        = $ARGV[2];
my $params_in        = $ARGV[3];
my $cpucount         = $ARGV[4];
my $PROCLU           = "$curdir/$ARGV[5]";
my $maxerror         = $ARGV[6];
my $maxflanconsid    = $ARGV[7];

# with dust filter (would be somewhat faster but less results)
#my $PROCLU_PARAM = "$curdir/eucledian.dst $curdir/eucledian.dst $params_in -p \"(0.70,0.30)\" -d -s $reffolder/reference.leb36";

# with PAM (doesnt take flanks into account) 
#my $PROCLU_PARAM = "$curdir/eucledian.dst $curdir/eucledian.dst $params_in -p \"(0.70,0.30)\" -s $reffolder/reference.leb36";

# basic
#my $PROCLU_PARAM = "$curdir/eucledian.dst $curdir/eucledian.dst $params_in -s $reffolder/reference.leb36";
my $PROCLU_PARAM_START = "$reffolder/reference.leb36";
my $PROCLU_PARAM_END   = "$curdir/eucledian.dst $params_in";

my $files_processed = 0;	# files processed
my %p; # associates forked pids with output pipe pids

my $MYLOCK=0;  # KA: meaningless, there is no shared memory in perl

# get a list of input files
opendir(DIR, $tgz_dir);
# the only extensions are .leb36
my @tarballs = grep(/\.(?:leb36)$/, readdir(DIR));
closedir(DIR);
my $tarball_count = @tarballs;
print "$tarball_count supported files found in $tgz_dir\n";
die "Exiting\n" if $tarball_count == 0;
$files_to_process = $tarball_count if $files_to_process > $tarball_count;

# enter dir
chdir($tgz_dir);


# fork as many new processes as there are CPUs
for (my $i = 0; $i < $cpucount; $i++) { $p{fork_proclu()} = 1 }

# wait for processes to finish and then fork new ones
while ((my $pid = wait) != -1) {

    # check return value
    my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
    if ($core) {
        warn "proclu process $pid dumped core\n";
        exit (1000);
    }
    elsif ($sig == 9) {
        warn "proclu process $pid was murdered!\n";
        exit (1001);
    }
    elsif ($rc != 0) {
        warn "proclu process $pid has returned $rc!\n";
        exit ($rc);
    }

    # KA: not sure why we're checking this
    die "Unrelated process $pid finished?\n" unless $p{$pid};

    # one instance has finished processing -- start a new one
    delete $p{$pid};
    $p{fork_proclu()} = 1;
}

print "Processing complete -- processed $files_processed file(s).\n";
warn strftime( "\nEnd: %F %T\n\n", localtime );

1;


############################ Procedures ###############################################################

sub fork_proclu {
    if ($files_processed >= $tarball_count) {
            return 0;
    }

    # this is meaningless, lock will just be memory duplicated
    # wait for shared variables to unlock
    while ($MYLOCK) { }

    # lock shared vars - no such thing
    $MYLOCK = 1;

    # use a predefined number of files
    my $until = $files_processed + $files_to_process - 1;
    $until = $tarball_count - 1 if $until > ($tarball_count - 1);
    my @file_slice = @tarballs[($files_processed)..($until)];
    my $file_slice_count = @file_slice;
    $files_processed += $file_slice_count;
    my $proclu_string;

    # unlock shared vars - still no such thing
    $MYLOCK = 0;

    defined(my $pid = fork)
        or die "Unable to fork: $!\n";

    # parent
    if ($pid != 0) { return $pid;}
    
    # child
    foreach (@file_slice) {
        $proclu_string = "$PROCLU $PROCLU_PARAM_START $_ $PROCLU_PARAM_END"
                         . "$maxerror $maxflanconsid 2>${_}.proclu_log";
        system($proclu_string);

        if ( $? == -1 ) { die "command failed: $!\n"; }
        else {
            my $rc = ($? >> 8);
            if ( 0 != $rc ) { print "proclu returned $rc ( $proclu_string  )!"; exit($rc); }
        }
    }

    # child must never return
    exit 0;
}


