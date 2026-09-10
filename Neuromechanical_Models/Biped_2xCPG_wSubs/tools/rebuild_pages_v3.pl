#!/usr/bin/perl
# Rebuild drawing pages v3 — single sequential pass, correct tails.
# For each page CDATA:
#   - collect element blocks (<Node>/<Link> drawings) with true end offsets
#   - nodes: emit verbatim, document order, dedup by Tag (last wins)
#   - links: dedup by Tag (last wins); recompute Org/Dst from the STANDARD link
#     endpoints (synapses AND adapter links) mapped to node ordinals on this page;
#     drop links whose endpoints aren't both drawn here or whose Tag has no standard link
#   - emit: prologue (with ShowGrid False + recomputed counts) + nodes + links + true tail
use strict; use warnings;
my ($in, $out) = @ARGV;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

# standard endpoints for ALL links (synapse + adapter)
my %end;    # link id -> [origin, dest]
while ($xml =~ m{<Link>(.*?)</Link>}gs) {
    my $b = $1;
    my ($id) = $b =~ /<ID>([^<]+)<\/ID>/;
    next unless $id;
    my ($o) = $b =~ /<OriginID>([^<]+)<\/OriginID>/;
    my ($d) = $b =~ /<DestinationID>([^<]+)<\/DestinationID>/;
    $end{$id} = [$o // '', $d // ''] if defined $o && defined $d;
}
printf "standard links with endpoints: %d\n", scalar(keys %end);

my $result = '';
my $pos = 0;
my ($totFixed, $totDropped, $totDedup) = (0, 0, 0);
while ($xml =~ /<DiagramXml><!\[CDATA\[/g) {
    my $cdata_start = pos($xml);
    my $cdata_end = index($xml, ']]></DiagramXml>', $cdata_start);
    last if $cdata_end < 0;
    my $d = substr($xml, $cdata_start, $cdata_end - $cdata_start);
    my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
    $pn =~ s/&amp;/&/g;

    # split: prologue (through "<AddFlow ...>"), body, epilogue (from "</AddFlow>")
    my $af = index($d, '<AddFlow ');
    die "no AddFlow in $pn\n" if $af < 0;
    my $body_start = index($d, '>', $af) + 1;
    my ($openTag) = $d =~ /(<AddFlow \b[^>]*>)/;
    my $epi = rindex($d, '</AddFlow>');
    my $prologue = substr($d, 0, $body_start);
    my $body     = substr($d, $body_start, $epi - $body_start);
    my $epilogue = substr($d, $epi);

    # walk element blocks with TRUE offsets
    my @nodes;   # [tag, text] in document order (dedup last-wins at emit)
    my @links;   # [tag, text] document order
    my ($lastEnd, $head) = (0, '');
    pos($body) = 0;
    my $guard = 0;
    while ($body =~ /<(Node|Link)\b/g) {
        my $etype = $1;
        my $estart = $-[0];
        $head = substr($body, 0, $estart) if !@nodes && !@links && $lastEnd == 0;
        my $closeTag = $etype eq 'Node' ? '</Node>' : '</Link>';
        my $eend = index($body, $closeTag, $estart);
        last if $eend < 0;
        $eend += length($closeTag);
        my $text = substr($body, $estart, $eend - $estart);
        my ($tag) = $text =~ /<Tag>([^<]*)<\/Tag>/;
        if ($etype eq 'Node') {
            my $t = defined($tag) ? $tag : "anon$guard";
            push @nodes, [$t, $text];
        } else {
            my $t = defined($tag) ? $tag : "anonL$guard";
            push @links, [$t, $text];
        }
        $guard++;
        pos($body) = $eend;
        $lastEnd = $eend;
    }
    my $trueTail = substr($body, $lastEnd);

    # dedup node tags (last wins) preserving first-occurrence position order? keep order of LAST
    my %nodeByText; my @nodeOrder;
    my %nodeSeen;
    for my $n (@nodes) {
        if ($nodeSeen{$n->[0]}++) {
            for my $e (@nodeOrder) { $e = $n if $e->[0] eq $n->[0]; }
        } else {
            push @nodeOrder, $n;
        }
    }
    # ordinals by EMITTED order (nodes first)
    my %ord;
    $ord{$nodeOrder[$_][0]} = $_ for 0 .. $#nodeOrder;

    # dedup link tags (last wins)
    my %linkLast;
    my @linkOrder;
    my %linkSeen;
    for my $l (@links) {
        if ($linkSeen{$l->[0]}++) {
            for my $e (@linkOrder) { $e = $l if $e->[0] eq $l->[0]; }
        } else {
            push @linkOrder, $l;
        }
    }

    # transform links
    my @emitLinks;
    my ($fixedHere, $droppedHere, $dedupHere) = (0, 0, 0);
    for my $l (@linkOrder) {
        my ($tag, $text) = @$l;
        my $isGen = $tag =~ /^b[01]a0/;
        unless (exists $end{$tag}) {
            # no standard endpoints known (junk drawing): keep originals ONLY if in range
            my ($org) = $text =~ /Org="(\d+)"/;
            my ($dst) = $text =~ /Dst="(\d+)"/;
            if (defined $org && defined $dst && $org < @nodeOrder && $dst < @nodeOrder) {
                push @emitLinks, $text;
            } else { $droppedHere++; }
            next;
        }
        my ($o, $dd) = @{$end{$tag}};
        if (exists $ord{$o} && exists $ord{$dd}) {
            my ($opentag) = $text =~ /^(<Link\b[^>]*>)/;
            my $rest = substr($text, length($opentag));
            $opentag =~ s/Org="\d+"/Org="$ord{$o}"/;
            $opentag =~ s/Dst="\d+"/Dst="$ord{$dd}"/;
            push @emitLinks, $opentag . $rest;
            $fixedHere++;
        } else {
            $droppedHere++;   # endpoints not drawn on this page
        }
    }
    $totFixed += $fixedHere; $totDropped += $droppedHere; $totDedup += $dedupHere;

    # prologue: grid off + exact counts
    my $newPro = $prologue;
    $newPro =~ s/<ShowGrid>True<\/ShowGrid>/<ShowGrid>False<\/ShowGrid>/;
    $newPro =~ s/(<AddFlow Nodes=")\d+(" Links=")\d+(")/$1 . scalar(@nodeOrder) . $2 . scalar(@emitLinks) . $3/e;

    $result .= substr($xml, $pos, $cdata_start - $pos);
    $result .= $newPro . $head . join('', map { $_->[1] } @nodeOrder) . join('', @emitLinks) . $trueTail . $epilogue;
    $pos = $cdata_end;
    print "page '$pn': nodes=".scalar(@nodeOrder)." links=".scalar(@emitLinks)." fixed=$fixedHere dropped=$droppedHere\n";
}
$result .= substr($xml, $pos);
print "SUMMARY: recomputed=$totFixed dropped=$totDropped\n";

open my $oh, '>', $out or die; print $oh $result; close $oh;
print "WROTE $out\n";
