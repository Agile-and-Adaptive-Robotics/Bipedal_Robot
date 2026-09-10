#!/usr/bin/perl
# Script A (v2): Complete RH-side wiring in Biped_2xCPG_wSubs.aproj
#  1. Canonically enumerate L-only synapses and mirror each to R (standard endpoints are truth)
#  2. Add commissural LH<->RH rhythm-generator connections (new OffPages on top page + 4 inhibitory synapses)
#  3. Give LH RG ext a 2 nA tonic kick (Deng)
#  4. Draw everything on the right CDATA pages
# Usage: perl build_A2.pl <in.aproj> <out.aproj>
use strict; use warnings;

my ($in, $out) = @ARGV;
die "usage: $0 in.aproj out.aproj\n" unless $in && $out;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

my $guid_seq = 0;
sub newguid { my $g = sprintf("b0a0%04d-0000-4000-8000-%012d", $guid_seq, 900000000000 + $guid_seq); $guid_seq++; return $g; }

# ---------- parse standard nodes (position scan, handles nesting truncation safely) ----------
my %node;
{
    pos($xml) = 0;
    while ($xml =~ /<Node>/g) {
        my $s = pos($xml);
        my $e = index($xml, '</Node>', $s);
        last if $e < 0;
        my $c = substr($xml, $s, $e - $s);
        my ($cls) = $c =~ /<ClassName>([^<]+)<\/ClassName>/;
        my ($id)  = $c =~ /<ID>([^<]+)<\/ID>/;
        my ($tx)  = $c =~ /<Text>([^<]*)<\/Text>/;
        if ($cls && $id) {
            my @p = split /\.|::/, $cls;
            my $short = $p[-1];
            if ($short =~ /^(NonSpiking|OffPage|Subsystem|StretchReceptor|Muscle|Spiking)$/) {
                my $lab = defined($tx) ? $tx : '';
                $lab =~ s/&amp;/&/g;
                $node{$id} = [$lab, $short];
            }
        }
        pos($xml) = $s;
    }
}
printf "parsed %d standard nodes\n", scalar keys %node;

# ---------- canonicalization (side-stripped role names) ----------
my %rolemap = (
  'Hip PF ext'=>'PFext','Hip PF flx'=>'PFflx','Hip PF ext IN'=>'INext','Hip PF flx IN'=>'INflx',
  'Hip MN ext'=>'MNext','Hip MN flx'=>'MNflx','Hip MN ext RE'=>'REext','Hip MN flx RE'=>'REflx',
  'Hip ext Ia'=>'SRiaExt','Hip flx Ia'=>'SRiaFlx','Hip ext Ia 2'=>'IaExt','Hip flx Ia 2'=>'IaFlx',
  'Hip ext Ib'=>'IbExt','Hip flx Ib'=>'IbFlx','HipZ Ext II'=>'IIext','HipZ FLX II'=>'IIflx',
  'Hip ext'=>'MUSext','Hip flx'=>'MUSflx',
  'HipZ PF E'=>'PFext','HipZ PF F'=>'PFflx','HipZ IN PF E'=>'INext','HipZ IN PF F'=>'INflx',
  'HipZ MN E'=>'MNext','HipZ MN F'=>'MNflx','HipZ R E'=>'REext','HipZ R F'=>'REflx',
  'HipZ-Ext Ia'=>'SRiaExt','HipZ-FLX Ia'=>'SRiaFlx','HipZ Ia E'=>'IaExt','HipZ Ia F'=>'IaFlx',
  'HipZ- Ext Ib'=>'IbExt','HipZ- Flx Ib'=>'IbFlx','HipZ-Ext II'=>'IIext','HipZ-FLX II'=>'IIflx',
  'Hip Ext'=>'MUSext','Hip Flx'=>'MUSflx',
  'RG ext'=>'RGext','RG flx'=>'RGflx','RG ext IN'=>'INRGext','RG flx IN'=>'INRGflx',
  'RG E'=>'RGext','RG F'=>'RGflx','IN RG E'=>'INRGext','IN RG F'=>'INRGflx',
);
sub canon {
    my $n = shift;
    (my $s = $n) =~ s/^(LH_|RH_|L |R )//;
    return $rolemap{$s} // $s;
}
sub sidelabel { my $n = shift; return ($n =~ /^(LH_|L )/) ? 'L' : (($n =~ /^(RH_|R )/) ? 'R' : 'X'); }

# ---------- mirror map: L endpoint -> R endpoint ----------
my %mirror = (
  '7635ff71'=>'095d4bd1-fe15-445b-80b8-459826254c8a',  # L RG ext real -> RH_RG E real
  'b82ffa13'=>'e04a7f4b-fcd5-4259-82f7-9ccc9d2747e2',  # L RG flx real -> RH_RG F real
  '2b22adea'=>'48485303-a756-4900-ae66-0d9881d68025',  # L RG ext IN -> RH_IN RG E
  'e211aaa9'=>'bbc3abd8-3f04-4343-b61e-e2d1c596abd2',  # L RG flx IN -> RH_IN RG F
  '116ae90c'=>'a422e14f-4f49-494f-a4de-c02cd6866cd1',  # L Hip PF ext real -> RH_HipZ PF E real
  'f519d63a'=>'dcff5f3e-f6a9-4351-92fc-c01f0b065a0b',  # L Hip PF flx real -> RH_HipZ PF F real
  '6088057a'=>'e10ffe67-c354-4c37-9075-bfbca1c132b9',  # L Hip PF ext IN -> RH_HipZ IN PF E
  'd9b2224c'=>'fb02aa52-e6ee-4cf3-a7ee-bf959ff433ea',  # L Hip PF flx IN -> RH_HipZ IN PF F
  '26d6abee'=>'88970959-393b-4221-87a8-3dea8e7381cc',  # L Hip MN ext -> RH_HipZ MN E
  'fb2595cc'=>'03573830-ba53-475e-b6de-efea78bf7c13',  # L Hip MN flx -> RH_HipZ MN F
  'a861b5cf'=>'956ae2d5-c2d5-4882-a1c9-994e7bb32e55',
  '7e2da9f7'=>'382fe4c3-7067-4cd7-a516-a9ddabb0199f',
  '396a9f7a'=>'755491d4-430a-42e9-86be-75b65678f015',
  '281f3bbb'=>'33680dd8-6aa2-458b-8f59-14887de4be1c',
  '256160de'=>'33680dd8-6aa2-458b-8f59-14887de4be1c',
  '1466cef2'=>'921a5273-4915-4ded-aa37-43972a8a6304',
  '927b20e0'=>'0bf31d1a-3932-42cc-9c04-048e092cd61e',
  '88d3b0f7'=>'7b7740cc-79fa-472e-9b49-37a1e630fc76',
  '8e5acd6f'=>'6b71d23d-7b9c-4280-9b0c-79bf8a7f519b',
  'e9b5192e'=>'9e2b2c4c-db4c-4dbb-a02f-c9b4f6cde9f2',
  'ca05607e'=>'62e62856-a392-4c32-afbf-69b98e3f4908',
  'dbb957f0'=>'a3b2a9a3-18a2-4f59-8e33-d9fb60a6440d',
  'cb722022'=>'29eafd95-9c6b-4d8c-b921-36865a89f7e4',
  'f7db5bec'=>'e5fcdd67-9e17-4309-81cd-bd8225dcb4f6',
  '8c219f3a'=>'466c7e5f-1c91-4e7b-8e65-641340d643a3',
  '3468b730'=>'4d637380-85aa-4176-ac77-389931eb7097',
  '6a7475fe'=>'c554929a-150b-4ede-ac4a-e71a5b3e9857',  # LH_AnkleZ- Ext Ib real -> RH real
  '04dc2dd1'=>'25fafd2b-d760-4fc1-9408-8fcbd8eec131',  # LH_HipZ Ext II real -> RH real
  '6abbc979'=>'5392c078-e469-46d6-99de-6f2ad335f07f',  # LH_HipZ FLX II real -> RH real
);

# ---------- enumerate synapses with canonical keys ----------
my @syn;    # {block, id, o, d, key, side}
{
    pos($xml) = 0;
    while ($xml =~ m{(<Link>.*?</Link>)}gs) {
        my $b = $1;
        next unless $b =~ /Behavior\.Synapse</;
        my ($id) = $b =~ /<ID>([^<]+)<\/ID>/;
        my ($o)  = $b =~ /<OriginID>([^<]+)<\/OriginID>/;
        my ($d)  = $b =~ /<DestinationID>([^<]+)<\/DestinationID>/;
        next unless $id && $o && $d && $node{$o} && $node{$d};
        my $so = sidelabel($node{$o}[0]);
        my $sd = sidelabel($node{$d}[0]);
        my $side = ($so eq 'L' && $sd ne 'R') ? 'L' : (($so eq 'R' && $sd ne 'L') ? 'R' : 'X');
        my $key = canon($node{$o}[0]).' => '.canon($node{$d}[0]);
        push @syn, { block=>$b, id=>$id, o=>$o, d=>$d, key=>$key, side=>$side };
    }
}
my (%Lkeys, %Rkeys);
for my $s (@syn) {
    if    ($s->{side} eq 'L') { $Lkeys{$s->{key}} = 1; }
    elsif ($s->{side} eq 'R') { $Rkeys{$s->{key}} = 1; }
}
my @missing = grep { !$Rkeys{$_} } sort keys %Lkeys;
print "canonical keys: L=".scalar(keys %Lkeys)." R=".scalar(keys %Rkeys)." L-only=".scalar(@missing)."\n";
print "  $_\n" for @missing;

# ---------- add_synapse machinery ----------
my @added_links;
sub add_to_nodelists {
    my ($xmlref, $nodeid, $linkid, $list) = @_;
    my $idx = index($$xmlref, "<ID>$nodeid</ID>");
    return unless $idx >= 0;
    my $nend = index($$xmlref, '</Node>', $idx);
    return if $nend < 0;
    my $seg_start = rindex($$xmlref, '<Node>', $idx);
    return if $seg_start < 0;
    my $seg = substr($$xmlref, $seg_start, $nend + 7 - $seg_start);
    my $newseg = $seg;
    if ($seg =~ /<$list>\s*<\/$list>/) { $newseg =~ s{<$list>\s*</$list>}{<$list>\n<ID>$linkid</ID>\n</$list>}; }
    elsif ($seg =~ /<\/$list>/)        { $newseg =~ s{</$list>}{<ID>$linkid</ID>\n</$list>}; }
    else { return; }
    substr($$xmlref, $seg_start, length($seg)) = $newseg;
}
sub add_synapse {
    my ($xmlref, $tmplblock, $src, $dst, $tag) = @_;
    my $nid = newguid();
    my $new = $tmplblock;
    $new =~ s/<ID>[^<]+<\/ID>/<ID>$nid<\/ID>/;
    $new =~ s{<OriginID>[^<]+</OriginID>}{<OriginID>$src</OriginID>} or die "no OriginID ($tag)";
    $new =~ s{<DestinationID>[^<]+</DestinationID>}{<DestinationID>$dst</DestinationID>} or die "no DestinationID ($tag)";
    add_to_nodelists($xmlref, $src, $nid, 'OutLinks');
    add_to_nodelists($xmlref, $dst, $nid, 'InLinks');
    my $pos = index($$xmlref, $tmplblock);
    die "template block not located\n" if $pos < 0;
    substr($$xmlref, $pos + length($tmplblock), 0) = "\n" . $new;
    print "ADDED $tag ($nid)\n";
    push @added_links, [$nid, $src, $dst, $tag];
}
sub r_link_exists {
    my ($o, $d) = @_;
    for my $s (@syn) {
        return 1 if $s->{side} eq 'R' && $s->{o} eq $o && $s->{d} eq $d;
    }
    return 0;
}

# ---------- STEP 1: mirror L-only canonical keys ----------
for my $key (@missing) {
    for my $s (@syn) {
        next unless $s->{side} eq 'L' && $s->{key} eq $key;
        my $no = $mirror{substr($s->{o},0,8)} or die "no mirror for origin $s->{o} ($key)\n";
        my $nd = $mirror{substr($s->{d},0,8)} or die "no mirror for dest $s->{d} ($key)\n";
        if (r_link_exists($no, $nd)) { print "SKIP (R has it): $key\n"; next; }
        add_synapse(\$xml, $s->{block}, $no, $nd, "MIRROR $key");
        my $nid = $added_links[-1][0];
        push @syn, { block=>'', id=>$nid, o=>$no, d=>$nd, key=>$key, side=>'R' };  # register
    }
}

# ---------- STEP 2: commissural LH<->RH RG connections ----------
my $offpage_tmpl;
{
    my $idx = index($xml, '<ID>256160de-7fd3-4594-a5a7-ee42c9a93e49</ID>');
    die "offpage template not found\n" if $idx < 0;
    my $s = rindex($xml, '<Node>', $idx);
    my $e = index($xml, '</Node>', $idx) + 7;
    $offpage_tmpl = substr($xml, $s, $e - $s);
    die "offpage template lacks Text\n" unless $offpage_tmpl =~ /<Text>/;
}
sub make_offpage {
    my ($linked, $text, $x, $y) = @_;
    my $n = $offpage_tmpl;
    my $id = newguid();
    $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
    $n =~ s/<Text>[^<]*<\/Text>/<Text>$text<\/Text>/;
    $n =~ s/<LinkedNodeID>[^<]*<\/LinkedNodeID>/<LinkedNodeID>$linked<\/LinkedNodeID>/;
    $n =~ s{<Location x="[^"]*" y="[^"]*"/>}{<Location x="$x" y="$y"/>};
    return ($id, $n);
}
my ($opLext, $b1) = make_offpage('7635ff71-07f1-485c-939a-42ba2f42714e', 'L RG ext',  80, 500);
my ($opLflx, $b2) = make_offpage('b82ffa13-18f8-4e68-a8f5-2c07edc4d8cf', 'L RG flx',  80, 580);
my ($opRext, $b3) = make_offpage('095d4bd1-fe15-445b-80b8-459826254c8a', 'RH_RG E', 700, 500);
my ($opRflx, $b4) = make_offpage('e04a7f4b-fcd5-4259-82f7-9ccc9d2747e2', 'RH_RG F', 700, 580);
{
    my $idx = index($xml, '<ID>9bf5105d-dfc9-4d76-83c4-542a765d1a7a</ID>');
    die "RH subsystem node not found\n" if $idx < 0;
    my $e = index($xml, '</Node>', $idx) + 7;
    substr($xml, $e, 0) = "\n$b1\n$b2\n$b3\n$b4";
    print "ADDED 4 commissural OffPage nodes (top level)\n";
}
my $comm_tmpl;
{ pos($xml) = 0;
  while ($xml =~ m{(<Link>.*?</Link>)}gs) { my $b=$1; if ($b =~ /Behavior\.Synapse</ && $b =~ /<ID>00ef6036/) { $comm_tmpl = $b; last; } }
  die "comm template not found\n" unless $comm_tmpl; }
add_synapse(\$xml, $comm_tmpl, $opLext, $opRext, 'COMMISSURAL L RG ext -| RH_RG E');
add_synapse(\$xml, $comm_tmpl, $opRext, $opLext, 'COMMISSURAL RH_RG E -| L RG ext');
add_synapse(\$xml, $comm_tmpl, $opLflx, $opRflx, 'COMMISSURAL L RG flx -| RH_RG F');
add_synapse(\$xml, $comm_tmpl, $opRflx, $opLflx, 'COMMISSURAL RH_RG F -| L RG flx');

# ---------- STEP 3: tonic kick on L RG ext ----------
{
    my $idx = index($xml, '<ID>7635ff71-07f1-485c-939a-42ba2f42714e</ID>');
    die "L RG ext not found\n" if $idx < 0;
    my $nend = index($xml, '</Node>', $idx);
    my $nstart = rindex($xml, '<Node>', $idx);
    my $seg = substr($xml, $nstart, $nend - $nstart);
    my $patched = $seg =~ s{<TonicStimulus Value="0" Scale="nano" Actual="0"/>}{<TonicStimulus Value="2" Scale="nano" Actual="2e-009"/>};
    die "tonic not patched in L RG ext\n" unless $patched == 1;
    substr($xml, $nstart, $nend - $nstart) = $seg;   # replace exactly the ORIGINAL range (seg grew)
    print "SET L RG ext TonicStimulus = 2 nA\n";
}

# ---------- STEP 4: CDATA drawings ----------
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
    my $old = $pages{$pn};
    my $pos = index($xml, $old);
    die "page $pn content not located\n" if $pos < 0;
    substr($xml, $pos, length($old)) = $newd;
    $pages{$pn} = $newd;
}
# draw the 4 new OffPage nodes on the top page (template: LH page's L RG ext OffPage drawing)
{
    my $lh = $pages{'LH'} or die "no LH page";
    my ($op_tmpl) = $lh =~ /(<Node\b[^>]*>(?:(?!<\/Node>).)*?<Tag>256160de[^<]*<\/Tag>(?:(?!<\/Node>).)*?<\/Node>)/s
        or die "OffPage drawing template not found on LH page";
    for my $o ([$opLext,'L RG ext',80,500], [$opLflx,'L RG flx',80,580], [$opRext,'RH_RG E',700,500], [$opRflx,'RH_RG F',700,580]) {
        my ($id,$txt,$x,$y) = @$o;
        my $n = $op_tmpl;
        $n =~ s{Left="[^"]*"}{Left="$x"}; $n =~ s{Top="[^"]*"}{Top="$y"};
        $n =~ s/<Text>[^<]*<\/Text>/<Text>$txt<\/Text>/s;
        $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$id<\/Tag>/s;
        my $d = $pages{'Neural Subsystem'};
        $d =~ s{<AddFlow Nodes="(\d+)" Links="(\d+)"}{'<AddFlow Nodes="'.($1+1).'" Links="'.$2.'"'}e;
        $d =~ s{(<\/AddFlow>)}{$n\n$1};
        replace_page('Neural Subsystem', $d);
    }
    print "DREW 4 OffPage nodes on top page\n";
}
my $global_link_tmpl;
{
    my $lh2 = $pages{'LH'};
    ($global_link_tmpl) = $lh2 =~ /(<Link\b[^>]*>(?:(?!<\/Link>).)*?<\/Link>)/s;
}
my $drawn = 0; my $undrawn = 0;
for my $al (@added_links) {
    my ($lid, $src, $dst, $tag) = @$al;
    my $host;
    for my $pn (keys %pages) {
        my $d = $pages{$pn};
        if (index($d, "<Tag>$src</Tag>") >= 0 && index($d, "<Tag>$dst</Tag>") >= 0) { $host = $pn; last; }
    }
    unless ($host) { $undrawn++; print "NOT DRAWN (functional only): $tag\n"; next; }
    my $d = $pages{$host};
    my ($tmpl) = $d =~ /(<Link\b[^>]*>(?:(?!<\/Link>).)*?<\/Link>)/s;
    $tmpl = $global_link_tmpl unless $tmpl;
    unless ($tmpl) { $undrawn++; print "no link template anywhere\n"; next; }
    my $newd = $tmpl;
    $newd =~ s{<Tag>[^<]*</Tag>}{<Tag>$lid</Tag>};
    $d =~ s{<AddFlow Nodes="(\d+)" Links="(\d+)"}{'<AddFlow Nodes="'.$1.'" Links="'.($2+1).'"'}e;
    $d =~ s{(<\/AddFlow>)}{$newd\n$1};
    replace_page($host, $d);
    $drawn++;
    print "DRAWN $tag on '$host'\n";
}

open my $ofh, '>', $out or die $!;
print $ofh $xml;
close $ofh;
print "WROTE $out\n";
print "SUMMARY: new synapses=".scalar(@added_links).", drawn=$drawn, functional-only=$undrawn\n";
