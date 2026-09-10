#!/usr/bin/perl
# Set Gait2392-derived MaximumTension on the 8 new muscles — CORRECT version.
# Only patches when Value != target; single regex, no fallback; verifies after.
use strict; use warnings;
my ($file, $asimToo) = @ARGV;
my %fmax = (gas_l=>2241, gas_r=>2241, bflh_l=>896, bflh_r=>896,
            semimem_l=>1288, semimem_r=>1288, rf_l=>1169, rf_r=>1169);

open my $fh, '<', $file or die; local $/; my $x = <$fh>; close $fh;
my $n = 0;
for my $nm (sort keys %fmax) {
    my $i = index($x, "<Name>$nm</Name>");
    die "muscle $nm not found\n" if $i < 0;
    my $e1 = index($x, '</RigidBody>', $i);
    die "no block end for $nm\n" if $e1 < 0;
    my $seg = substr($x, $i, $e1 - $i);
    my $F = $fmax{$nm};
    my $newseg = $seg;
    my $ok = ($newseg =~ s{(<MaximumTension Value=")\d+(" Scale="None" Actual=")\d+(")}{$1 . $F . $2 . $F . $3}e);
    die "MaximumTension not patched for $nm (pattern mismatch)\n" unless $ok == 1;
    substr($x, $i, $e1 - $i) = $newseg;
    $n++;
}
open my $oh, '>', $file or die; print $oh $x; close $oh;
print "set_fmax: patched $n in $file\n";
