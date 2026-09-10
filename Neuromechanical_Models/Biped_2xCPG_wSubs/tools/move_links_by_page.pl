#!/usr/bin/perl
# Move generated links (b0a0/b1a0) into the fragment whose page draws them,
# and turn off ShowGrid on all diagram pages.
use strict; use warnings;
my ($in, $out) = @ARGV;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

sub find_close {
    my ($strref, $start, $tag) = @_;
    my $close = "</$tag>";
    my $openExact = "<$tag>";
    my $openAttr  = "<$tag ";
    my $depth = 1;
    my $p = $start;
    while (1) {
        my $nc = index($$strref, $close, $p);
        return -1 if $nc < 0;
        my $no1 = index($$strref, $openExact, $p);
        my $no2 = index($$strref, $openAttr, $p);
        my $no = -1;
        $no = $no1 if $no1 >= 0;
        $no = $no2 if $no2 >= 0 && ($no < 0 || $no2 < $no);
        # CDATA content is opaque text: skip it entirely
        my $cd = index($$strref, '<![CDATA[', $p);
        if ($cd >= 0 && $cd < $nc && ($no < 0 || $cd < $no)) {
            my $ce = index($$strref, ']]>', $cd);
            return -1 if $ce < 0;
            $p = $ce + 3;
            next;
        }
        if ($no >= 0 && $no < $nc) { $depth++; $p = $no + length($openExact); next; }
        $depth--;
        return $nc + length($close) if $depth == 0;
        $p = $nc + length($close);
    }
}

# ---- enumerate Subsystem fragments (span + Text) ----
my @frags;   # {start, end, id, name}
{
    my $p = 0;
    while ((my $i = index($xml, 'Behavior.Nodes.Subsystem', $p)) >= 0) {
        my $open = rindex($xml, '<Node>', $i);
        if ($open < 0) { $p = $i + 10; next; }
        # skip if this span already recorded (child subsystems share the parent scan region)
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
print "fragments:\n";
for my $fr (sort { $a->{start} <=> $b->{start} } @frags) {
    print "  [$fr->{start}..$fr->{end}] $fr->{name} ($fr->{id})\n";
}

sub innermost_frag {
    my ($pos) = @_;
    my $best;
    for my $fr (@frags) {
        next unless $pos > $fr->{start} && $pos < $fr->{end};
        $best = $fr if !defined($best) || $fr->{end} - $fr->{start} < $best->{end} - $best->{start};
    }
    return $best;
}

# ---- map page name -> fragment ----
my %pageFrag;
for my $fr (@frags) { $pageFrag{$fr->{name}} = $fr if $fr->{name}; }

# ---- for each generated standard link: find its drawing page; move if needed ----
my @genIds;
while ($xml =~ /<ID>(b[01]a0[^<]+)<\/ID>\s*<AssemblyFile>IntegrateFireGUI\.dll<\/AssemblyFile>\s*<ClassName>IntegrateFireGUI\.DataObjects\.Behavior\.Synapse</g) {
    push @genIds, $1;
}
# fallback: any standard synapse Link whose ID starts with b0a0/b1a0
{
    my $p2 = 0; my %seen;
    while ($xml =~ m{<Link>(.*?)</Link>}gs) {
        my $b = $1;
        if ($b =~ /Behavior\.Synapse</ && $b =~ /<ID>(b[01]a0[^<]+)<\/ID>/ && !$seen{$1}++) { push @genIds, $1; }
    }
}
my @genIdsU = do { my %s; grep { !$s{$_}++ } @genIds };
print "generated standard synapse links: " . scalar(@genIdsU) . "\n";

my $moved = 0; my $already = 0;
for my $lid (sort @genIdsU) {
    # find drawing page: the CDATA containing <Tag>lid</Tag>
    my $dpos = index($xml, "<Tag>$lid</Tag>");
    if ($dpos < 0) { print "  $lid: NOT DRAWN anywhere (functional only)\n"; next; }
    # which page CDATA contains dpos
    my $pageName;
    for my $fr2 (0) { }
    {
        # find the nearest preceding PageName whose CDATA start precedes dpos
        pos($xml) = 0;
        while ($xml =~ /<DiagramXml><!\[CDATA\[/g) {
            my $s = pos($xml);
            my $e = index($xml, ']]></DiagramXml>', $s);
            last if $e < 0;
            if ($dpos > $s && $dpos < $e) {
                my $d = substr($xml, $s, $e - $s);
                my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
                $pn =~ s/&amp;/&/g;
                $pageName = $pn;
                last;
            }
            pos($xml) = $s;
        }
    }
    die "no page for drawing of $lid\n" unless $pageName;
    my $tfrag = $pageFrag{$pageName} or die "no fragment for page $pageName\n";
    # current block position + enclosing fragment
    my $i = index($xml, $tfrag ? "<ID>$lid</ID>" : "", 0);
    # find the standard link block: <Link>...<ID>lid...
    my $bpos;
    {
        pos($xml) = 0;
        while ($xml =~ /<Link>/g) {
            my $s = pos($xml);
            my $e = index($xml, '</Link>', $s);
            last if $e < 0;
            my $c = substr($xml, $s, $e - $s);
            if ($c =~ /<ID>\Q$lid\E<\/ID>/) { $bpos = $s - 6; last; }
            pos($xml) = $s;
        }
    }
    die "standard block for $lid not found\n" unless defined $bpos;
    my $cur = innermost_frag($bpos);
    my $curName = $cur ? $cur->{name} : '(top)';
    if ($cur && $cur->{name} eq $pageName) { $already++; next; }
    # move: remove block, insert before target fragment close
    my $blk = substr($xml, $bpos, index($xml, '</Link>', $bpos) + 7 - $bpos);
    substr($xml, $bpos, length($blk)) = '';
    # target fragment close may have shifted; recompute by fragment node id
    my $tclose = do {
        my $idi = index($xml, "<ID>$tfrag->{id}</ID>");
        my $open = rindex($xml, '<Node>', $idi);
        find_close(\$xml, $open + 6, 'Node') - 7;
    };
    substr($xml, $tclose, 0) = "\n" . $blk;
    $moved++;
    print "MOVED $lid: $curName -> $pageName\n";
}
print "moved=$moved already-ok=$already\n";

# ---- grid off on all pages ----
my $g = ($xml =~ s/<ShowGrid>True<\/ShowGrid>/<ShowGrid>False<\/ShowGrid>/g);
print "grids turned off: $g\n";

open my $oh, '>', $out or die; print $oh $xml; close $oh;
print "WROTE $out\n";
