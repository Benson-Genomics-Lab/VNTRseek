#!/usr/bin/env perl

use strict;
use warnings;

use FileHandle;
use Getopt::Std;
use File::Copy;

my %options;
getopts( 'rn', \%options );

# like the shell getopt, "d:" means d takes an argument
die "Usage:
	./renumber.pl [-r -n] [directory]

Renumbers pairs of LEB36 and index files.

\t-r\tno renumbering -- only sign is switched
\t-n\tno indexes -- only work on LEB36 files
" if not scalar(@ARGV);

my $switch_sign = ( defined $options{'r'} && $options{'r'} ) ? 1 : 0;
my $use_index   = ( defined $options{'n'} && $options{'n'} ) ? 0 : 1;

print STDERR "Only changing the sign\n" if $switch_sign;
print STDERR "Not using indexes\n"      if !$use_index;

my @files;

if ( -d $ARGV[0] ) {

    # is a directory
    opendir( DIR, $ARGV[0] ) or die "Cannot open $ARGV[0] directory\n";
    my @dir_content = readdir(DIR);
    closedir(DIR);

    my @set_leb36 = sort grep {s/\.leb36$//} @dir_content;

    if ($use_index) {
        my @set_index = sort grep {s/\.index$//} @dir_content;

        # now find symmetric difference between @set_leb36 and @set_index
        my %count;
        my @symm_diff;
        foreach ( @set_leb36, @set_index ) { $count{$_}++ }
        foreach ( keys %count ) {
            if ( $count{$_} < 2 ) {
                push @symm_diff, $_;
            }
        }
        die 'No valid pairs for ' . join( ', ', @symm_diff ) . "\n"
            if scalar(@symm_diff);
    }

    @files = map { $ARGV[0] . '/' . $_ } @set_leb36;
    my $file_count = scalar(@files);
    print STDERR "$file_count LEB36 files found in directory $ARGV[0]\n";
}

my $UNR = 0;    # "universal number of records

my $output_index;
my $output_index_h;
my $input_index_h;
my $index_line;
my @index_line;

if ($switch_sign) {
    foreach my $input (@files) {
        my $output_leb36   = "$input.leb36.tmp";
        my $output_leb36_h = FileHandle->new(">$output_leb36")
            or die "Cannot open $output_leb36 for writing\n";
        my $input_leb36_h = FileHandle->new("<$input.leb36")
            or die "Cannot open $input.leb36 for reading\n";

        if ($use_index) {
            $output_index   = "$input.index.tmp";
            $output_index_h = FileHandle->new(">$output_index")
                or die "Cannot open $output_index for writing\n";
            $input_index_h = FileHandle->new("<$input.index")
                or die "Cannot open $input.index for reading\n";
        }

        while ( my $leb36_line = <$input_leb36_h> ) {
            $UNR++;
            my @leb36_line = split /\s+/, $leb36_line;  # split on white space

            if ($use_index) {
                $index_line = <$input_index_h>;
                chomp $index_line;
                die "Too few lines in $input.index\n"
                    if not defined $index_line;
                @index_line = split /\t+/, $index_line; # split on white space
                die
                    "Numbers differ at line $UNR in $input.index and $input.leb36\n"
                    if $leb36_line[0] != $index_line[0];
                $index_line[0] = -abs( $index_line[0] );
                print $output_index_h join( "\t", @index_line ) . "\n";
            }

            $leb36_line[0] = -abs( $leb36_line[0] );
            print $output_leb36_h join( ' ', @leb36_line ) . "\n";
        }
        if ($use_index) {
            die "Too many lines in $input.index\n"
                if defined <$input_index_h>;
            $output_index_h->close
                or die "Cannot close $output_index file handle\n";
            $input_index_h->close
                or die "Cannot close $input.index file handle\n";
            move( $output_index, "$input.index" );
        }

        $output_leb36_h->close
            or die "Cannot close $output_leb36 file handle\n";
        $input_leb36_h->close
            or die "Cannot close $input.leb36 file handle\n";
        move( $output_leb36, "$input.leb36" );
    }
}
else {
    foreach my $input (@files) {
        my $output_leb36   = "$input.leb36.tmp";
        my $output_leb36_h = FileHandle->new(">$output_leb36")
            or die "Cannot open $output_leb36 for writing\n";
        my $input_leb36_h = FileHandle->new("<$input.leb36")
            or die "Cannot open $input.leb36 for reading\n";

        if ($use_index) {
            $output_index   = "$input.index.tmp";
            $output_index_h = FileHandle->new(">$output_index")
                or die "Cannot open $output_index for writing\n";
            $input_index_h = FileHandle->new("<$input.index")
                or die "Cannot open $input.index for reading\n";
        }

        while ( my $leb36_line = <$input_leb36_h> ) {
            $UNR++;
            if ( $UNR <= 0 ) {
                $output_leb36_h->close;
                $input_leb36_h->close;
                unlink $output_leb36;
                if ($use_index) {
                    $output_index_h->close;
                    $input_index_h->close;
                    unlink $output_index;
                }
                die "Integer overflow\n";
            }
            my @leb36_line = split /\s+/, $leb36_line;  # split on white space

            if ($use_index) {
                $index_line = <$input_index_h>;
                chomp $index_line;
                die "Too few lines in $input.index\n"
                    if not defined $index_line;
                @index_line = split /\t+/, $index_line; # split on white space
                die
                    "Numbers differ at line $UNR in $input.index and $input.leb36\n"
                    if $leb36_line[0] != $index_line[0];
                $index_line[0] = $UNR;
                print $output_index_h join( "\t", @index_line ) . "\n";
            }

            $leb36_line[0] = $UNR;
            print $output_leb36_h join( ' ', @leb36_line ) . "\n";
        }
        if ($use_index) {
            die "Too many lines in $input.index\n"
                if defined <$input_index_h>;
            $output_index_h->close
                or die "Cannot close $output_index file handle\n";
            $input_index_h->close
                or die "Cannot close $input.index file handle\n";
            move( $output_index, "$input.index.renumbered" );
        }
        $output_leb36_h->close
            or die "Cannot close $output_leb36 file handle\n";
        $input_leb36_h->close
            or die "Cannot close $input.leb36 file handle\n";
        move( $output_leb36, "$input.leb36.renumbered" );
    }
}

print STDERR "$UNR records processed\n";
0;
