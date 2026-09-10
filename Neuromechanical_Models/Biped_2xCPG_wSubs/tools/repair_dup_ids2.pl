#!/usr/bin/perl
# Correct repair v2: re-GUID duplicated child-object IDs (occurrences 2..n) using a
# single global substitution with /ge — no positional bookkeeping to get wrong.
use strict; use warnings;
my ($in, $out) = @ARGV;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

my $seq = 0;
sub ng { return sprintf("b3a0%04d-0000-4000-8000-%012d", $seq, 700000000000 + $seq++); }
my $CHILD = qr(<(?:StimulusTension|LengthTension|Gain|CaActivation|CaDeactivation|PID)>\s*<ID>);

# pass 1: count
my %count;
while ($xml =~ /$CHILD([^<]+)<\/ID>/g) { $count{$1}++; }

# pass 2: rename 2nd+ occurrences
my %seen; my $fixed = 0;
$xml =~ s{($CHILD)([^<]+)(</ID>)}{
    my ($pre, $id, $post) = ($1, $2, $3);
    if ($count{$id} > 1) {
        $seen{$id}++;
        if ($seen{$id} > 1) { $fixed++; $pre . ng() . $post; }
        else { $pre . $id . $post; }
    } else {
        $pre . $id . $post;
    }
}ge;

open my $ofh, '>', $out or die $!;
print $ofh $xml;
close $ofh;
print "renamed $fixed duplicated child IDs\n";
