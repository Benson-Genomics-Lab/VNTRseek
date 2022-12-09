package Product;

use strict;
use warnings;
use v5.24;
use List::Util 'all';

sub new {
    my ($class, $reps, @items) = @_;
    my @index = (0) x $reps;
    $index[-1] = -1;

    my $self = bless {
        reps => $reps,
        index => \@index,
        items => \@items,
        last => scalar @items - 1
    }, $class;
}

sub next {
    my $self = shift;
    if (all { $_ == $self->{last}} @{$self->{index}}) { return; }

    for (my $i = $self->{reps} - 1; $i >= 0; $i--) {
        if (@{$self->{index}}[$i] == @{$self->{items}} - 1) {
            @{$self->{index}}[$i] = 0
        }
        else {
            @{$self->{index}}[$i]++;
            last;
        }
    }
    return @{$self->{items}}[@{$self->{index}}];
}

1;