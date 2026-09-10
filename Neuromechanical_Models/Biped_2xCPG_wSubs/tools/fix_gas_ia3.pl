#!/usr/bin/perl
# Surgical fix v3: Gas IaE -> antagonist-Ia synapse targeted the HIP Ia F instead of the
# ANKLE Ia F. Repoint DestinationID on both sides, move InLinks bookkeeping, and
# re-add the page drawings that fix_orgdst removed for those links.
use strict; use warnings;
my ($in, $out) = @ARGV;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

sub find_link_block {   # standard <Link> block whose own <ID> starts with prefix
    my $p = shift;
    pos($xml) = 0;
    while ($xml =~ m{(<Link>.*?</Link>)}gs) {
        my $b = $1;
        my ($id) = $b =~ /<ID>([^<]+)<\/ID>/;
        next unless defined $id;
        return $b if $id =~ /^\Q$p\E/;
    }
    die "link block '$p' not found\n";
}
sub find_node_block {   # <Node> block whose own <ID> equals the full id
    my $id = shift;
    pos($xml) = 0;
    while ($xml =~ /<Node>/g) {
        my $s = pos($xml);
        my $e = index($xml, '</Node>', $s);
        last if $e < 0;
        my $c = substr($xml, $s, $e - $s);
        my ($m) = $c =~ /<ID>([^<]+)<\/ID>/;
        if (defined $m && $m eq $id) {
            my $before = substr($xml, 0, $s);
            my $start = rindex($before, '<Node>') >= 0 ? rindex($before, '<Node>') : $s;
            return (substr($xml, $start, $e + 7 - $start), $start, $e + 7);
        }
        pos($xml) = $s;
    }
    die "node $id not found\n";
}
sub fullid {   # full GUID for an 8-char prefix, from node or link definitions
    my $p = shift;
    pos($xml) = 0;
    while ($xml =~ /<Node>/g) {
        my $s = pos($xml);
        my $e = index($xml, '</Node>', $s);
        last if $e < 0;
        my $c = substr($xml, $s, $e - $s);
        my ($m) = $c =~ /<ID>([^<]+)<\/ID>/;
        if (defined $m && $m =~ /^\Q$p\E/) { pos($xml) = $s; return $m; }
        pos($xml) = $s;
    }
    my $b = find_link_block($p);
    my ($id) = $b =~ /<ID>([^<]+)<\/ID>/;
    return $id;
}

my @fixes = (
  ['b1a00051', 'baecf9d6', 'LH_Anklez Motoneuron'],   # L_Gas Ia E -> LH_AnkleZ Ia F
  ['b1a00149', '11ed4906', 'RH_Anklez Motoneuron'],   # RH_Gas Ia E -> RH_AnkleZ Ia F
);
my @redraw;   # [linkFullId, page, link8]
for my $f (@fixes) {
    my ($lid8, $newDest8, $page) = @$f;
    my $lid    = fullid($lid8);
    my $newDst = fullid($newDest8);
    my $blk = find_link_block($lid8);
    my ($oldDst) = $blk =~ /<DestinationID>([^<]+)<\/DestinationID>/;
    die "no destination in $lid8\n" unless $oldDst;
    my $nblk = $blk;
    $nblk =~ s{<DestinationID>\Q$oldDst\E</DestinationID>}{<DestinationID>$newDst</DestinationID>} or die "dest repoint failed";
    my $pos = index($xml, $blk);
    die "link block not located\n" if $pos < 0;
    substr($xml, $pos, length($blk)) = $nblk;
    # bookkeeping: drop stale InLinks entry on old dest; add on new dest
    my ($oseg, $os, $oe) = find_node_block($oldDst);
    my $newoseg = $oseg;
    $newoseg =~ s{<ID>\Q$lid\E</ID>\s*}{};
    substr($xml, $os, $oe - $os) = $newoseg if $newoseg ne $oseg;
    my ($nseg, $ns, $ne) = find_node_block($newDst);
    my $newnseg = $nseg;
    if ($newnseg =~ /<InLinks>\s*<\/InLinks>/) {
        $newnseg =~ s{<InLinks>\s*</InLinks>}{<InLinks>\n<ID>$lid</ID>\n</InLinks>};
    } elsif ($newnseg =~ /<\/InLinks>/) {
        $newnseg =~ s{</InLinks>}{<ID>$lid</ID>\n</InLinks>};
    }
    substr($xml, $ns, $ne - $ns) = $newnseg;
    print "REPOINTED $lid8: dest repointed\n";
    push @redraw, [$lid, $page, $lid8];
}

# redraw on the ankle pages with node-ordinal Org/Dst
my %pages;
{
    my $work = $xml;
    while ($work =~ /<DiagramXml><!\[CDATA\[/g) {
        my $s = pos($work);
        my $e = index($work, ']]></DiagramXml>', $s);
        last if $e < 0;
        my $d = substr($work, $s, $e - $s);
        my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
        $pn =~ s/&amp;/&/g;
        $pages{$pn} = $d;
        pos($work) = $s;
    }
}
sub replace_page {
    my ($pn, $newd) = @_;
    my $pos = index($xml, $pages{$pn});
    die "page $pn not located\n" if $pos < 0;
    substr($xml, $pos, length($pages{$pn})) = $newd;
    $pages{$pn} = $newd;
}
for my $r (@redraw) {
    my ($lid, $pn, $lid8) = @$r;
    my $d = $pages{$pn};
    my $blk = find_link_block($lid8);
    my ($o)   = $blk =~ /<OriginID>([^<]+)</;
    my ($dst) = $blk =~ /<DestinationID>([^<]+)</;
    my (%ord, $k);
    $k = 0; pos($d) = 0;
    while ($d =~ /<Node\b/g) {
        my $ns = pos($d); my $ne = index($d, '</Node>', $ns); last if $ne < 0;
        my $c = substr($d, $ns, $ne - $ns);
        my ($t) = $c =~ /<Tag>([^<]*)<\/Tag>/;
        $ord{$t} = $k++ if $t;
        pos($d) = $ns;
    }
    die "endpoint $o not drawn on $pn\n" unless exists $ord{$o};
    die "endpoint $dst not drawn on $pn\n" unless exists $ord{$dst};
    my ($tmpl) = $d =~ /(<Link\b[^>]*>(?:(?!<\/Link>).)*?<\/Link>)/s;
    my $n = $tmpl;
    $n =~ s/Org="\d+"/Org="$ord{$o}"/;
    $n =~ s/Dst="\d+"/Dst="$ord{$dst}"/;
    $n =~ s{<Tag>[^<]*</Tag>}{<Tag>$lid</Tag>};
    $d =~ s{<AddFlow Nodes="(\d+)" Links="(\d+)"}{'<AddFlow Nodes="'.$1.'" Links="'.($2+1).'"'}e;
    $d =~ s{(<\/AddFlow>)}{$n\n$1};
    replace_page($pn, $d);
    print "REDREW $lid8 on '$pn' (org=$ord{$o} dst=$ord{$dst})\n";
}

open my $ofh, '>', $out or die $!;
print $ofh $xml;
close $ofh;
print "WROTE $out\n";
