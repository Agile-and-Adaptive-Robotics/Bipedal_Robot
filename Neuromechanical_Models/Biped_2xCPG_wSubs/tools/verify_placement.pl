#!/usr/bin/perl
# Verify: every object drawn on a page lives in that page's fragment.
use strict; use warnings;
my $file = shift;
open my $fh, '<', $file or die $!; local $/; my $xml = <$fh>; close $fh;

sub find_close {
    my ($strref, $start, $tag) = @_;
    my $close = "</$tag>";
    my $depth = 1;
    my $p = $start;
    while (1) {
        my $no = index($$strref, "<$tag>", $p);
        my $nc = index($$strref, $close, $p);
        return -1 if $nc < 0;
        if ($no >= 0 && $no < $nc) { $depth++; $p = $no + length("<$tag>"); next; }
        $depth--;
        return $nc + length($close) if $depth == 0;
        $p = $nc + length($close);
    }
}

# fragments
my @frags;
{
    my $p = 0;
    while ((my $i = index($xml, 'Behavior.Nodes.Subsystem', $p)) >= 0) {
        my $open = rindex($xml, '<Node>', $i);
        my $e = find_close(\$xml, $open + 6, 'Node');
        my ($id) = substr($xml, $open, $e - $open) =~ /<ID>([^<]+)<\/ID>/;
        my ($nm) = substr($xml, $open, $e - $open) =~ /<Text>([^<]*)<\/Text>/;
        $nm = defined($nm) ? $nm : '';
        $nm =~ s/&amp;/&/g;
        push @frags, { start => $open, end => $e, id => $id, name => $nm }
            unless grep { $_->{start} == $open } @frags;
        $p = $i + 10;
    }
}
sub innermost_name {
    my ($pos) = @_;
    my $best;
    for my $fr (@frags) {
        next unless $pos > $fr->{start} && $pos < $fr->{end};
        $best = $fr if !defined($best) || $fr->{end} - $fr->{start} < $best->{end} - $best->{start};
    }
    return $best ? $best->{name} : '(none)';
}

# object id -> its fragment name + type
my %obj;    # id -> {frag, kind}
pos($xml) = 0;
while ($xml =~ /<Node>/g) {
    my $s = pos($xml);
    my $e = find_close(\$xml, $s, 'Node');
    last if $e < 0;
    my $blk = substr($xml, $s - 6, $e - $s + 6);
    my ($id) = $blk =~ /<ID>([^<]+)<\/ID>/;
    if ($id) {
        my $frag = innermost_name($s);
        my ($cls) = $blk =~ /<ClassName>([^<]+)<\/ClassName>/;
        my @pp = split /\.|::/, $cls // '';
        $obj{$id} = { frag => $frag, kind => $pp[-1] };
    }
    pos($xml) = $s;
}

# pages: drawn objects vs their home fragments
my $issues = 0;
pos($xml) = 0;
while ($xml =~ /<DiagramXml><!\[CDATA\[/g) {
    my $s = pos($xml);
    my $e = index($xml, ']]></DiagramXml>', $s);
    last if $e < 0;
    my $d = substr($xml, $s, $e - $s);
    my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
    $pn =~ s/&amp;/&/g;
    my (%nodeTags, @linkDrawn);
    while ($d =~ /<Node\b[^>]*>(.*?)<\/Node>/gs) {
        my ($t) = $1 =~ /<Tag>([^<]+)<\/Tag>/;
        $nodeTags{$t} = 1 if $t;
    }
    while ($d =~ /<Link\b[^>]*>(.*?)<\/Link>/gs) {
        my ($t) = $1 =~ /<Tag>([^<]+)<\/Tag>/;
        push @linkDrawn, $t if $t;
    }
    my ($nodeIssues, $linkIssues) = (0, 0);
    for my $t (sort keys %nodeTags) {
        next unless $obj{$t};
        if ($obj{$t}{frag} ne $pn) { $nodeIssues++; print "  [$pn] node '$t' lives in '$obj{$t}{frag}'\n" if $nodeIssues <= 3; }
    }
    for my $t (@linkDrawn) {
        next unless $obj{$t};
        if ($obj{$t}{frag} ne $pn) { $linkIssues++; print "  [$pn] link '$t' lives in '$obj{$t}{frag}'\n" if $linkIssues <= 3; }
    }
    printf "%-28s nodes=%d links=%d nodeIssues=%d linkIssues=%d\n", $pn, scalar(keys %nodeTags), scalar(@linkDrawn), $nodeIssues, $linkIssues;
    $issues += $nodeIssues + $linkIssues;
}
print "TOTAL placement issues: $issues\n";
