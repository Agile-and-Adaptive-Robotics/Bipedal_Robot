#!/usr/bin/perl
# Script B (v2): Add gastrocnemius, biceps femoris long head, rectus femoris, semimembranosus
# to each leg of Biped_2xCPG_wSubs.aproj (physical + neural, Deng-style wiring).
# Usage: perl build_B2.pl <in.aproj> <out.aproj>
use strict; use warnings;

my ($in, $out) = @ARGV;
die "usage\n" unless $in && $out;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

my $gseq = 0;
sub ng { return sprintf("b1a0%04d-0000-4000-8000-%012d", $gseq, 800000000000 + $gseq++); }

# ---------------- generic helpers ----------------
my %node;
{
    pos($xml) = 0;
    while ($xml =~ /<Node>/g) {
        my $s = pos($xml);
        my $e = index($xml, '</Node>', $s);
        last if $e < 0;
        my $c = substr($xml, $s, $e - $s);
        my ($id) = $c =~ /<ID>([^<]+)<\/ID>/;
        my ($tx) = $c =~ /<Text>([^<]*)<\/Text>/;
        $node{$id} = defined($tx) ? $tx : '' if $id;
        pos($xml) = $s;
    }
}
sub find_syn {
    my ($o8, $d8) = @_;
    pos($xml) = 0;
    while ($xml =~ m{(<Link>.*?</Link>)}gs) {
        my $b = $1;
        next unless $b =~ /Behavior\.Synapse</;
        my ($o) = $b =~ /<OriginID>([^<]+)<\/OriginID>/;
        my ($d) = $b =~ /<DestinationID>([^<]+)<\/DestinationID>/;
        return $b if $o =~ /^\Q$o8\E/ && $d =~ /^\Q$d8\E/;
    }
    die "synapse $o8->$d8 not found\n";
}
sub find_adlink {
    my ($o8, $d8) = @_;
    pos($xml) = 0;
    while ($xml =~ m{(<Link>.*?</Link>)}gs) {
        my $b = $1;
        next unless $b =~ /Behavior\.Links\.Adapter</;
        my ($o) = $b =~ /<OriginID>([^<]+)<\/OriginID>/;
        my ($d) = $b =~ /<DestinationID>([^<]+)<\/DestinationID>/;
        return $b if $o =~ /^\Q$o8\E/ && $d =~ /^\Q$d8\E/;
    }
    die "adapter link $o8->$d8 not found\n";
}
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
sub insert_after_block {
    my ($xmlref, $block, $new) = @_;
    my $pos = index($$xmlref, $block);
    die "block not located for insert\n" if $pos < 0;
    substr($$xmlref, $pos + length($block), 0) = "\n" . $new;
}
sub resolve_full {   # expand 8-char GUID prefix to full id
    my $id = shift;
    return $id if length($id) > 8;
    my $i = index($xml, "<ID>$id");
    die "cannot resolve prefix $id\n" if $i < 0;
    my $s = $i + 4; my $e = index($xml, '</ID>', $s);
    return substr($xml, $s, $e - $s);
}
sub clone_link {
    my ($tmplblock, $o, $d) = @_;
    $o = resolve_full($o); $d = resolve_full($d);
    my $nid = ng();
    my $new = $tmplblock;
    $new =~ s/<ID>[^<]+<\/ID>/<ID>$nid<\/ID>/;
    $new =~ s{<OriginID>[^<]+</OriginID>}{<OriginID>$o</OriginID>} or die "no OriginID";
    $new =~ s{<DestinationID>[^<]+</DestinationID>}{<DestinationID>$d</DestinationID>} or die "no DestinationID";
    add_to_nodelists(\$xml, $o, $nid, 'OutLinks');
    add_to_nodelists(\$xml, $d, $nid, 'InLinks');
    insert_after_block(\$xml, $tmplblock, $new);
    return $nid;
}

# ---------------- drawing helpers ----------------
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
    die "page $pn not located\n" if $pos < 0;
    substr($xml, $pos, length($old)) = $newd;
    $pages{$pn} = $newd;
}
my $drk = 0;
sub draw_node {
    my ($pn, $id, $text) = @_;
    my $d = $pages{$pn} or die "no page $pn";
    my ($tmpl) = $d =~ /(<Node\b[^>]*>(?:(?!<\/Node>).)*?<\/Node>)/s;
    die "no node drawing template on $pn\n" unless $tmpl;
    my $n = $tmpl;
    my ($tl) = $n =~ /Left="([-\d.]+)"/; my ($tt) = $n =~ /Top="([-\d.]+)"/;
    my $nx = ($tl // 100) + 35 + ($drk % 3) * 20;
    my $ny = ($tt // 100) + 28 + int($drk / 3) * 26; $drk++;
    $n =~ s{Left="[^"]*"}{Left="$nx"}; $n =~ s{Top="[^"]*"}{Top="$ny"};
    $n =~ s/<Text>[^<]*<\/Text>/<Text>$text<\/Text>/s;
    $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$id<\/Tag>/s;
    $d =~ s{<AddFlow Nodes="(\d+)" Links="(\d+)"}{'<AddFlow Nodes="'.($1+1).'" Links="'.$2.'"'}e;
    $d =~ s{(<\/AddFlow>)}{$n\n$1};
    replace_page($pn, $d);
}
my $global_link_tmpl;
{ my $lh = $pages{'LH'}; ($global_link_tmpl) = $lh =~ /(<Link\b[^>]*>(?:(?!<\/Link>).)*?<\/Link>)/s; }
sub draw_link {
    my ($pn, $lid) = @_;
    my $d = $pages{$pn} or die "no page $pn";
    my ($tmpl) = $d =~ /(<Link\b[^>]*>(?:(?!<\/Link>).)*?<\/Link>)/s;
    $tmpl = $global_link_tmpl unless $tmpl;
    my $n = $tmpl;
    $n =~ s{<Tag>[^<]*</Tag>}{<Tag>$lid</Tag>};
    $d =~ s{<AddFlow Nodes="(\d+)" Links="(\d+)"}{'<AddFlow Nodes="'.$1.'" Links="'.($2+1).'"'}e;
    $d =~ s{(<\/AddFlow>)}{$n\n$1};
    replace_page($pn, $d);
}

# ---------------- config ----------------
my %cfg = (
 L => {
  P=>'L',
  page=>{ankle=>'LH_Anklez Motoneuron', hip=>'LH_HipZ', knee=>'LH_Knee Motoneuron'},
  mnE_ankle=>'fbf80f0e', mnF_ankle=>'3468b730', muscleNode_ankle=>'beeb9a51', srNode_ankle=>'23f2ce3c',
  n2p_ankle=>'70437f90', p2nIa_ankle=>'aa043851', p2nIb_ankle=>'f40920d0',
  affIa_ankle=>'9e2d3758', iaE_ankle=>'985125fc', iaF_ankle=>'baecf9d6', ib_ankle=>'6a7475fe',
  mnE_hip=>'26d6abee', mnF_hip=>'fb2595cc', muscleNode_hip=>'bc198dd5', srNode_hip=>'17b7cb57',
  n2p_hip=>'d81b03a8', p2nIa_hip=>'a21975cb', p2nIb_hip=>'7086c1b0',
  affIa_hip=>'0a7518d5', iaE_hip=>'3a622fba', iaF_hip=>'3ee81205', ib_hip=>'1abd3af9',
  pfE_hip_real=>'116ae90c', pfF_hip_real=>'f519d63a',
  kaE_real=>'9ceb5a3d', kaF_real=>'5f4c2b51', kaF_offpage_anklepg=>"dfc28c45", kaE_offpage_anklepg=>"548c27d1",
  kneeMNE=>'e7f0f45b', kneeMNF=>'8b3d2c0c',
  offpageTmpl=>'256160de',
  phys=>{
    thigh_muscle=>'c7d73eb3', thigh_rec=>'473396c9', calf_muscle=>'9765f0b1', calf_rec=>'7157203f',
    att_thighO=>'cb3c439f', att_thighI=>'7a174733', att_calfI=>'720c245f',
    att_hipO_ext=>'a056d4c3', att_hipO_flex=>'399da34b',
  },
 },
 R => {
  P=>'RH',
  page=>{ankle=>'RH_Anklez Motoneuron', hip=>'RH_HipZ', knee=>'RH_Knee Motoneuron'},
  mnE_ankle=>'5f4e432e', mnF_ankle=>'4d637380', muscleNode_ankle=>'95c11a83', srNode_ankle=>'2b13ffc8',
  n2p_ankle=>'7d633c8b', p2nIa_ankle=>'811e7a10', p2nIb_ankle=>'5955ca33',
  affIa_ankle=>'c65326b4', iaE_ankle=>'9a7c7f49', iaF_ankle=>'11ed4906', ib_ankle=>'c554929a',
  mnE_hip=>'88970959', mnF_hip=>'03573830', muscleNode_hip=>'89bae014', srNode_hip=>'fba15b80',
  n2p_hip=>'904b97eb', p2nIa_hip=>'c3d700cd', p2nIb_hip=>'653f6f6e',
  affIa_hip=>'d4c4ff77', iaE_hip=>'be4136ad', iaF_hip=>'ce567ed6', ib_hip=>'77da980e',
  pfE_hip_real=>'a422e14f', pfF_hip_real=>'dcff5f3e',
  kaE_real=>'b86372f3', kaF_real=>'58a2e50e', kaF_offpage_anklepg=>"e19f57dd", kaE_offpage_anklepg=>"f7eed3a5",
  kneeMNE=>'22970258', kneeMNF=>'716e9049',
  offpageTmpl=>'256160de',
  phys=>{
    thigh_muscle=>'ab982010', thigh_rec=>'dda95abe', calf_muscle=>'6014988f', calf_rec=>'eef3765e',
    att_thighO=>'?knee_R_flx_origin', att_thighI=>'?knee_R_flx_insertion', att_calfI=>'?ankle_R_ext_insertion',
    att_hipO_ext=>'?hip_R_ext_origin', att_hipO_flex=>'?hip_R_flx_origin',
  },
 },
);
# resolve ?name attachment ids
my %attIdByName;
{
    pos($xml) = 0;
    while ($xml =~ /<RigidBody>/g) {
        my $s = pos($xml);
        my $e = index($xml, '</RigidBody>', $s);
        last if $e < 0;
        my $c = substr($xml, $s, $e - $s);
        my ($nm) = $c =~ /<Name>([^<]*)<\/Name>/;
        my ($id) = $c =~ /<ID>([^<]+)<\/ID>/;
        my ($tp) = $c =~ /<Type>([^<]*)<\/Type>/;
        $attIdByName{$nm} = $id if $nm && $id && ($tp||'') eq 'Attachment';
        pos($xml) = $s;
    }
}
for my $S (values %cfg) {
    for my $k (keys %{$S->{phys}}) {
        my $v = $S->{phys}{$k};
        if ($v =~ /^\?(.+)$/) { my $id = $attIdByName{$1} or die "attachment '$1' not found\n"; $S->{phys}{$k} = $id; }
    }
}

# ---------------- muscle specs ----------------
# role: 'ext' (extends hip/ankle: Gas,BFlh,Semimem) or 'flx' (RF at hip)
my @muscles = (
  { key=>'Gas',     ss=>'ankle', role=>'ext', maxT=>1500, rest=>44, drive=>"KAE",   sameAg=>'mnF_ankle', iaXJoint=>"kneeMNF" },
  { key=>"BFlh",    ss=>"hip",   role=>"ext", maxT=>1500, rest=>34, drive=>"hipE",  sameAg=>"mnF_hip",   iaXJoint=>"kneeMNF" },
  { key=>"Semimem", ss=>"hip",   role=>"ext", maxT=>1500, rest=>34, drive=>"hipE",  sameAg=>"mnF_hip",   iaXJoint=>"kneeMNF" },
  { key=>'RF',      ss=>'hip',   role=>'flx', maxT=>1500, rest=>34, drive=>'hipF',  sameAg=>'mnE_hip',   iaXJoint=>'kneeMNF' },
);

# ================= STEP 0: fix pre-existing R-side Renshaw-drive =================
# RH knee/ankle N2P muscle-drive adapters were wired from Renshaw cells (RE/RF) instead of MNs.
# Repoint them to the corresponding MNs so RH mirrors LH (hip) and the paper's architecture.
{
    my @fix = (
      ['23bd37e2', '42f67e06', '22970258'],   # R knee ext: R E  -> RH_Knee MN E
      ['b4e2d77e', '1c45a020', '716e9049'],   # R knee flx: R F  -> RH_Knee MN F
      ['7d633c8b', '5a1fec02', '5f4e432e'],   # R ankle ext: R E -> RH_AnkleZ MN E
      ['5fb72a05', '42778e81', '4d637380'],   # R ankle flx: R F -> RH_AnkleZ MN F
    );
    for my $f (@fix) {
        my ($n2p, $oldO, $newO) = @$f;
        my $full = sub { my $p = shift; my $i = index($xml, "<ID>$p"); die "id $p not found\n" if $i < 0; my $s2 = $i + 4; my $e2 = index($xml, '</ID>', $s2); return substr($xml, $s2, $e2 - $s2); };
        my ($fN2p, $fOld, $fNew) = ($full->($n2p), $full->($oldO), $full->($newO));
        # 1) the adapter node's OriginID
        my $idx = index($xml, "<ID>$fN2p</ID>");
        my $s = rindex($xml, '<Node>', $idx);
        my $e = index($xml, '</Node>', $idx) + 7;
        my $seg = substr($xml, $s, $e - $s);
        my $nseg = $seg;
        $nseg =~ s{<OriginID>\Q$fOld\E</OriginID>}{<OriginID>$fNew</OriginID>} or die "N2P $n2p origin not repointed";
        substr($xml, $s, $e - $s) = $nseg;
        # 2) the adapter link feeding it
        pos($xml) = 0;
        my $done = 0;
        while ($xml =~ m{(<Link>.*?</Link>)}gs) {
            my $b = $1;
            next unless $b =~ /Behavior\.Links\.Adapter</;
            my ($o) = $b =~ /<OriginID>([^<]+)<\/OriginID>/;
            my ($d) = $b =~ /<DestinationID>([^<]+)<\/DestinationID>/;
            next unless (defined $o && $o eq $fOld) && (defined $d && $d eq $fN2p);
            my ($lid) = $b =~ /<ID>([^<]+)<\/ID>/;
            my $nb = $b;
            $nb =~ s{<OriginID>\Q$fOld\E</OriginID>}{<OriginID>$fNew</OriginID>};
            substr($xml, index($xml, $b), length($b)) = $nb;
            add_to_nodelists(\$xml, $fNew, $lid, 'OutLinks');
            print "REPOINTED Renshaw-drive adapter $n2p: origin -> $newO\n";
            $done = 1;
            last;
        }
        die "adapter link into $n2p not found\n" unless $done;
    }
}

# ================= STEP 1: physical =================
my %physId;
for my $side ('L','R') {
    my $S = $cfg{$side};
    my $lp = lc $side;
    for my $m (@muscles) {
        my $k = $m->{key};
        my $calf = ($k eq 'Gas') ? 1 : 0;
        # attachments (positions copied exactly from templates; no offsets)
        my ($oTmpl, $oName, $oParent, $iTmpl, $iName, $iParent);
        if ($calf) {
            ($oTmpl,$oName,$oParent) = ($S->{phys}{att_thighO}, "gas_${lp}_origin", ($side eq 'L' ? 'femur_L' : 'femur_R'));
            ($iTmpl,$iName,$iParent) = ($S->{phys}{att_calfI},  "gas_${lp}_insertion", ($side eq 'L' ? 'foot_L' : 'foot_R'));
        } elsif ($k eq 'RF') {
            ($oTmpl,$oName,$oParent) = ($S->{phys}{att_hipO_flex}, "rf_${lp}_origin", 'ROOT');
            ($iTmpl,$iName,$iParent) = ($S->{phys}{att_thighI}, "rf_${lp}_insertion", ($side eq 'L' ? 'tibia_L' : 'tibia_R'));
        } else {
            my $base = ($k eq 'BFlh') ? "bflh_${lp}" : "semimem_${lp}";
            ($oTmpl,$oName,$oParent) = ($S->{phys}{att_hipO_ext}, "${base}_origin", 'ROOT');
            ($iTmpl,$iName,$iParent) = ($S->{phys}{att_thighI}, "${base}_insertion", ($side eq 'L' ? 'tibia_L' : 'tibia_R'));
        }
        my %att;
        for my $a ([$oTmpl,$oName,$oParent,'origin'], [$iTmpl,$iName,$iParent,'insertion']) {
            my ($tid,$tname,$parent,$role) = @$a;
            my $idx = index($xml, "<ID>$tid");
            die "attachment template $tid not found\n" if $idx < 0;
            my $s2 = rindex($xml, '<RigidBody>', $idx);
            my $e2 = index($xml, '</RigidBody>', $idx) + 12;
            my $n = substr($xml, $s2, $e2 - $s2);
            my $id = ng();
            $n =~ s/<Name>[^<]*<\/Name>/<Name>$tname<\/Name>/;
            $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
            my $cbIdx;
            if ($parent eq 'ROOT') {
                my $ref = index($xml, '<Name>hip_L_flx_rec</Name>');
                $ref = index($xml, '<Name>hip_R_flx_rec</Name>') if $ref < 0;
                $cbIdx = rindex($xml, '<ChildBodies>', $ref);
            } else {
                my $ref = index($xml, "<Name>$parent</Name>");
                die "parent $parent not found\n" if $ref < 0;
                $cbIdx = index($xml, '<ChildBodies>', $ref);
            }
            die "no ChildBodies for $parent\n" if $cbIdx < 0;
            substr($xml, $cbIdx + 13, 0) = "\n" . $n;
            $att{$role} = $id;
        }
        # muscle + receptor bodies (into ROOT ChildBodies)
        my $physNm = lc($k) . "_${lp}";
        my $mTmpl = $calf ? $S->{phys}{calf_muscle} : $S->{phys}{thigh_muscle};
        my $rTmpl = $calf ? $S->{phys}{calf_rec}    : $S->{phys}{thigh_rec};
        for my $spec ([$mTmpl, $physNm, $m->{maxT}, $m->{rest}, 'muscle'], [$rTmpl, "${physNm}_rec", 100, $m->{rest}, 'rec']) {
            my ($tid,$tname,$tmax,$trest,$kind) = @$spec;
            my $idx = index($xml, "<ID>$tid");
            my $s2 = rindex($xml, '<RigidBody>', $idx);
            my $e2 = index($xml, '</RigidBody>', $idx) + 12;
            my $n = substr($xml, $s2, $e2 - $s2);
            my $id = ng();
            $n =~ s/<Name>[^<]*<\/Name>/<Name>$tname<\/Name>/;
            $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
            my $seen = 0;
            $n =~ s{<AttachID>[^<]*</AttachID>}{ $seen++ == 0 ? "<AttachID>$att{origin}</AttachID>" : "<AttachID>$att{insertion}</AttachID>" }ge;
            die "template $tid had ".(0+$seen)." attachids (need 2)\n" unless $seen == 2;
            if ($kind eq 'muscle') {
                my $c1 = $n =~ s{(<MaximumTension Value=")[^"]*(" Scale="None" Actual=")[^"]*(")}{$1$tmax$2$tmax$3};
                die "MaximumTension not patched for $tname\n" unless $c1 == 1;
            }
            my $c2 = $n =~ s{(<RestingLength Value=")[^"]*(" Scale="centi" Actual=")[^"]*(")}{$1 . $trest . $2 . ($trest*0.01) . $3}e;
            die "RestingLength not patched for $tname (scale not centi?)\n" unless $c2 == 1;
            my $ref = index($xml, '<Name>hip_L_flx_rec</Name>');
            $ref = index($xml, '<Name>hip_R_flx_rec</Name>') if $ref < 0;
            my $cbIdx = rindex($xml, '<ChildBodies>', $ref);
            substr($xml, $cbIdx + 13, 0) = "\n" . $n;
            $physId{"$side.$k"}{$kind} = $id;
        }
        $physId{"$side.$k"}{origin} = $att{origin};
        $physId{"$side.$k"}{insertion} = $att{insertion};
        print "PHYS $side $k done (muscle $physId{\"$side.$k\"}{muscle})\n";
    }
}

# ================= STEP 2: neural =================
my (@newNodes, @newLinks);
sub clone_neuron_node {
    my ($tmplId, $newName) = @_;
    my $idx = index($xml, "<ID>$tmplId");
    my $s = rindex($xml, '<Node>', $idx);
    my $e = index($xml, '</Node>', $idx) + 7;
    my $blk = substr($xml, $s, $e - $s);
    my $id = ng();
    my $n = $blk;
    $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
    $n =~ s/<Text>[^<]*<\/Text>/<Text>$newName<\/Text>/;
    $n =~ s{<InLinks>.*?</InLinks>}{<InLinks/>}gs;
    $n =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks/>}gs;
    insert_after_block(\$xml, $blk, $n);
    return $id;
}
sub clone_part_node {   # Muscle / StretchReceptor neural node
    my ($tmplId, $newName, $physPartId) = @_;
    my $idx = index($xml, "<ID>$tmplId");
    my $s = rindex($xml, '<Node>', $idx);
    my $e = index($xml, '</Node>', $idx) + 7;
    my $blk = substr($xml, $s, $e - $s);
    my $id = ng();
    my $n = $blk;
    $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
    $n =~ s/<Text>[^<]*<\/Text>/<Text>$newName<\/Text>/;
    $n =~ s/<LinkedBodyPartID>[^<]*<\/LinkedBodyPartID>/<LinkedBodyPartID>$physPartId<\/LinkedBodyPartID>/;
    $n =~ s{<InLinks>.*?</InLinks>}{<InLinks/>}gs;
    $n =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks/>}gs;
    insert_after_block(\$xml, $blk, $n);
    return $id;
}
sub clone_adapter_node {
    my ($tmplId, $newName, $o, $d) = @_;
    my $idx = index($xml, "<ID>$tmplId");
    my $s = rindex($xml, '<Node>', $idx);
    my $e = index($xml, '</Node>', $idx) + 7;
    my $blk = substr($xml, $s, $e - $s);
    my $id = ng();
    my $n = $blk;
    $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
    $n =~ s/<Text>[^<]*<\/Text>/<Text>$newName<\/Text>/;
    $n =~ s{<OriginID>[^<]+</OriginID>}{<OriginID>$o</OriginID>} or die "adapter lacks OriginID";
    $n =~ s{<DestinationID>[^<]+</DestinationID>}{<DestinationID>$d</DestinationID>} or die "adapter lacks DestinationID";
    $n =~ s{<InLinks>.*?</InLinks>}{<InLinks/>}gs;
    $n =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks/>}gs;
    insert_after_block(\$xml, $blk, $n);
    return $id;
}
sub clone_offpage {
    my ($tmplId, $linkedId, $text, $insertAfterId) = @_;
    $linkedId = resolve_full($linkedId);
    my $idx = index($xml, "<ID>$tmplId");
    my $s = rindex($xml, '<Node>', $idx);
    my $e = index($xml, '</Node>', $idx) + 7;
    my $blk = substr($xml, $s, $e - $s);
    my $id = ng();
    my $n = $blk;
    $n =~ s/<ID>[^<]+<\/ID>/<ID>$id<\/ID>/;
    $n =~ s/<Text>[^<]*<\/Text>/<Text>$text<\/Text>/;
    $n =~ s/<LinkedNodeID>[^<]*<\/LinkedNodeID>/<LinkedNodeID>$linkedId<\/LinkedNodeID>/;
    $n =~ s{<InLinks>.*?</InLinks>}{<InLinks/>}gs;
    $n =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks/>}gs;
    # insert after a node of the TARGET subsystem so page ownership is right
    my $aIdx = index($xml, "<ID>$insertAfterId");
    my $aEnd = index($xml, '</Node>', $aIdx) + 7;
    substr($xml, $aEnd, 0) = "\n" . $n;
    return $id;
}

my %neural;
for my $side ('L','R') {
    my $S = $cfg{$side};
    my $P = $S->{P};
    my $adv = 41;
    for my $m (@muscles) {
        my $k = $m->{key};
        my $ankle = ($m->{ss} eq 'ankle') ? 1 : 0;
        my $pg  = $ankle ? $S->{page}{ankle} : $S->{page}{hip};
        my $ph  = $physId{"$side.$k"};
        my $mn  = clone_neuron_node($ankle ? $S->{mnF_ankle} : $S->{mnE_hip}, "${P}_${k} MN");
        my $aff = clone_neuron_node($ankle ? $S->{affIa_ankle} : $S->{affIa_hip}, "${P}_${k}-Ext Ia");
        my $iaE = clone_neuron_node($ankle ? $S->{iaE_ankle} : $S->{iaE_hip}, "${P}_${k} Ia E");
        my $ib  = clone_neuron_node($ankle ? $S->{ib_ankle} : $S->{ib_hip}, "${P}_${k}- Ext Ib");
        my $mus = clone_part_node($ankle ? $S->{muscleNode_ankle} : $S->{muscleNode_hip}, "${P}_${k}", $ph->{muscle});
        my $sr  = clone_part_node($ankle ? $S->{srNode_ankle} : $S->{srNode_hip}, "${P}_${k} ", $ph->{rec});
        my $adv0 = $adv;
        my $a1 = clone_adapter_node($ankle ? $S->{n2p_ankle} : $S->{n2p_hip}, $adv0,   $mn, $mus);
        my $a2 = clone_adapter_node($ankle ? $S->{p2nIa_ankle} : $S->{p2nIa_hip}, $adv0+1, $sr, $aff);
        my $a3 = clone_adapter_node($ankle ? $S->{p2nIb_ankle} : $S->{p2nIb_hip}, $adv0+2, $mus, $ib);
        $adv += 3;
        $neural{"$side.$k"} = { mn=>$mn, aff=>$aff, iaE=>$iaE, ib=>$ib, mus=>$mus, sr=>$sr };
        push @newNodes, map { [$pg, @$_] } (
          [$mn,"${P}_${k} MN"],[$aff,"${P}_${k}-Ext Ia"],[$iaE,"${P}_${k} Ia E"],[$ib,"${P}_${k}- Ext Ib"],
          [$mus,"${P}_${k}"],[$sr,"${P}_${k} "],[$a1,"$adv0"],[$a2,($adv0+1).""],[$a3,($adv0+2).""]);
        # adapter links
        my ($mET, $n2T, $muT, $srT, $iaT, $affT, $ibT, $libT);
        if ($ankle) { ($mET,$n2T,$muT,$srT,$iaT,$affT,$ibT,$libT) =
          ($S->{mnE_ankle},$S->{n2p_ankle},$S->{muscleNode_ankle},$S->{srNode_ankle},$S->{p2nIa_ankle},$S->{affIa_ankle},$S->{p2nIb_ankle},$S->{ib_ankle}); }
        else        { ($mET,$n2T,$muT,$srT,$iaT,$affT,$ibT,$libT) =
          ($S->{mnE_hip},$S->{n2p_hip},$S->{muscleNode_hip},$S->{srNode_hip},$S->{p2nIa_hip},$S->{affIa_hip},$S->{p2nIb_hip},$S->{ib_hip}); }
        my $id1 = clone_link(find_adlink($mET, $n2T), $mn,  $a1);
        my $id2 = clone_link(find_adlink($n2T, $muT), $a1,  $mus);
        my $id3 = clone_link(find_adlink($srT, $iaT), $sr,  $a2);
        my $id4 = clone_link(find_adlink($iaT, $affT), $a2, $aff);
        my $id5 = clone_link(find_adlink($muT, $ibT), $mus, $a3);
        my $id6 = clone_link(find_adlink($ibT, $libT), $a3, $ib);
        push @newLinks, map { [$pg, $_] } ($id1,$id2,$id3,$id4,$id5,$id6);

        # ------- synapses -------
        # drive
        if ($ankle) {
            my $srcOff = ($m->{drive} eq "KAE") ? $S->{kaE_offpage_anklepg} : $S->{kaF_offpage_anklepg};
            my $dstTmpl = ($m->{drive} eq "KAE") ? $S->{mnE_ankle} : $S->{mnF_ankle};
            my $lid = clone_link(find_syn($srcOff, $dstTmpl), $srcOff, $mn);
            push @newLinks, [$pg, $lid];
        } else {
            my ($srcPF, $dstMNt);
            if ($m->{drive} eq 'hipF') { ($srcPF,$dstMNt) = ($S->{pfF_hip_real}, $S->{mnF_hip}); }
            else                       { ($srcPF,$dstMNt) = ($S->{pfE_hip_real}, $S->{mnE_hip}); }
            my $lid = clone_link(find_syn($srcPF, $dstMNt), $srcPF, $mn);
            push @newLinks, [$pg, $lid];
            my $kaReal = ($m->{drive} eq 'hipF') ? $S->{kaE_real} : $S->{kaF_real};
            my $opText = ($m->{drive} eq "hipF") ? "${P}_K&amp;A PF E" : "${P}_K&amp;A PF F";
            my $opid = clone_offpage($S->{offpageTmpl}, $kaReal, $opText, $S->{mnE_hip});
            push @newNodes, [$pg, $opid, $opText];
            my $lid2 = clone_link(find_syn($S->{kaF_offpage_anklepg}, $S->{mnF_ankle}), $opid, $mn);
            push @newLinks, [$pg, $lid2];
        }
        # afferent -> IaE
        {
            my ($o,$d) = $ankle ? ($S->{affIa_ankle},$S->{iaE_ankle}) : ($S->{affIa_hip},$S->{iaE_hip});
            my $lid = clone_link(find_syn($o,$d), $aff, $iaE);
            push @newLinks, [$pg, $lid];
        }
        # IaE -> same-subsystem antagonist MN
        {
            my $agMn = $S->{$m->{sameAg}};
            my ($o,$d) = $ankle ? ($S->{iaE_ankle},$S->{mnF_ankle}) : ($S->{iaE_hip},$S->{mnF_hip});
            my $lid = clone_link(find_syn($o,$d), $iaE, $agMn);
            push @newLinks, [$pg, $lid];
        }
        # IaE -> antagonist-side Ia proper (same subsystem)
        my ($s4o, $s4d);
        if ($ankle)        { ($s4o,$s4d) = ($S->{iaE_ankle}, $S->{iaF_ankle}); }
        elsif ($m->{role} eq 'ext') { ($s4o,$s4d) = ($S->{iaE_hip}, $S->{iaF_hip}); }
        else                        { ($s4o,$s4d) = ($S->{iaF_hip}, $S->{iaE_hip}); }   # RF: flexor-side Ia inhibits extensor Ia
        {
            my $tgt = ($m->{role} eq 'flx') ? $S->{iaE_hip} : $S->{iaF_hip};
            my $lid = clone_link(find_syn($s4o,$s4d), $iaE, $tgt);
            push @newLinks, [$pg, $lid];
        }
        # antagonist-side Ia -> new MN (disynaptic reciprocal)
        {
            my ($o,$d);
            if ($ankle)                { ($o,$d) = ($S->{iaF_ankle}, $S->{mnE_ankle}); }   # flexor Ia -> extensor MN
            elsif ($m->{role} eq 'ext'){ ($o,$d) = ($S->{iaF_hip}, $S->{mnE_hip}); }
            else                       { ($o,$d) = ($S->{iaE_hip}, $S->{mnF_hip}); }       # RF: extensor Ia -> flexor MN
            my $lid = clone_link(find_syn($o,$d), $o, $mn);
            push @newLinks, [$pg, $lid];
        }
        # cross-joint Ia via OffPage on knee page
        {
            my $kneePg = $S->{page}{knee};
            my $kneeMN = $S->{$m->{iaXJoint}};
            my $opid = clone_offpage($S->{offpageTmpl}, $iaE, "${P}_${k} Ia E", $S->{kneeMNF});
            push @newNodes, [$kneePg, $opid, "${P}_${k} Ia E"];
            my ($o,$d) = $ankle ? ($S->{iaE_ankle},$S->{mnF_ankle}) : ($S->{iaE_hip},$S->{mnF_hip});
            my $lid = clone_link(find_syn($o,$d), $opid, $kneeMN);
            push @newLinks, [$kneePg, $lid];
        }
        # Ib -> own MN
        {
            my ($o,$d) = $ankle ? ($S->{ib_ankle},$S->{mnE_ankle}) : ($S->{ib_hip},$S->{mnE_hip});
            my $lid = clone_link(find_syn($o,$d), $ib, $mn);
            push @newLinks, [$pg, $lid];
        }
    }
}

# cleanup: strip stale In/OutLinks on script-A OffPages (b0a0 ids)
{
    pos($xml) = 0;
    my $fixed = 0;
    while ($xml =~ /<Node>/g) {
        my $s = pos($xml);
        my $e = index($xml, '</Node>', $s);
        last if $e < 0;
        my $c = substr($xml, $s, $e - $s);
        if ($c =~ /<ID>b0a0/ && $c =~ /<LinkedNodeID>/) {   # script-A OffPage
            my $n = $c;
            $n =~ s{<InLinks>.*?</InLinks>}{<InLinks/>}gs;
            $n =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks/>}gs;
            if ($n ne $c) { substr($xml, $s, $e - $s) = $n; $fixed++; }
        }
        pos($xml) = $s;
    }
    print "cleaned $fixed script-A OffPages\n";
}

# ================= STEP 3: drawings =================
for my $n (@newNodes) { draw_node(@$n); }
for my $l (@newLinks) { draw_link(@$l); }

open my $ofh, '>', $out or die $!;
print $ofh $xml;
close $ofh;
printf "WROTE $out : newNodes=%d newLinks=%d\n", scalar(@newNodes), scalar(@newLinks);
