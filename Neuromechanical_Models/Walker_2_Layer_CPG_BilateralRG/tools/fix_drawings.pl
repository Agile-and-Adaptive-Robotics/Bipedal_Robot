#!/usr/bin/perl
# fix_drawings.pl — rebuild ALL appended page drawings in
# Walker_2_Layer_CPG_BilateralRG.aproj with correct AddFlow endpoints,
# proper shapes, and deliberate node placement.
#
# Verified convention (verify_handles.pl, 2026-09-16): a drawing <Link>'s
# Org/Dst = 0-based index of the endpoint NODE entry in the CDATA's
# interleaved (nodes+links) file order. All 147 original synapse drawings
# match exactly; the 62 appended ones kept template endpoints (invisible
# duplicates) — the bug Ben saw.
#
# Operations:
#   1. delete my 78 appended link drawings (Tag starts 'cafe')
#   2. delete the 4 orphan crossed-link drawings (functional links removed
#      in build_rg; drawings kept back then because deletion seemed to break
#      handles — safe now since ALL later drawings get recomputed)
#   3. retemplate my 8 adapter node drawings to the real adapter look
#   4. add 4 node drawings for the contact rigid bodies (their incoming
#      adapter links need a drawn endpoint; originals draw bodies too)
#   5. place all 40 new node drawings in deliberate, collision-checked
#      groups (R RG mirrors L RG; c1/V3 beside it; II near Ia relays;
#      contact row at bottom)
#   6. re-add all 78 link drawings with computed Org/Dst, no stale Points
#   7. fix AddFlow counters, validate
use strict; use warnings;
my $dir = shift @ARGV or die "usage: fix_drawings.pl <dir>\n";
my $PROJ = "$dir/Walker_2_Layer_CPG_BilateralRG.aproj";
sub slurp { open my $f,'<',$_[0] or die "$_[0]: $!"; local $/; my $s=<$f>; close $f; $s }
sub spew  { open my $o,'>',$_[0] or die "$_[0]: $!"; print {$o} $_[1]; close $o }
sub g { my $n=shift; sprintf("cafe%04d-0000-4000-8000-%012d",$n,900000000000+$n) }
sub resolve { my ($x,$p)=@_; my $i=index($x,"<ID>$p"); die "prefix $p not found\n" if $i<0;
  substr($x,$i+4,index($x,'</ID>',$i+4)-$i-4) }

my $x = slurp($PROJ);
my $CS = index($x,'<![CDATA[')+9;
my $CE = index($x,']]></DiagramXml>',$CS);
my $cd = substr($x,$CS,$CE-$CS);

# ---------- endpoint map: functional link id -> [originId, destId] ----------
my %func;
pos($x)=0;
while ($x =~ /<Link>(.*?)<\/Link>/gs) {
  my $b=$1;
  next unless $b=~/Behavior\.Synapse</ || $b=~/Behavior\.Links\.Adapter</;
  my ($id)=$b=~/<ID>([^<]+)<\/ID>/; my ($o)=$b=~/<OriginID>([^<]+)<\/OriginID>/; my ($d)=$b=~/<DestinationID>([^<]+)<\/DestinationID>/;
  next unless $id && $o && $d;
  $func{$id}=[$o,$d];
}

# ---------- parse CDATA into entry list ----------
# entry: {kind=>'N'|'L', tag, text, block, start(len in cd)} — we work on a
# list and rebuild the CDATA at the end.
my @ents;
pos($cd)=0;
while ($cd =~ /(<Node\s+Left=[^>]*>|<Link\s+Org=[^>]*>)/g) {
  my $open=$1; my $kind=$open=~/^<Node/?"N":"L";
  my $from=pos($cd);
  my $close=$kind eq "N"?"</Node>":"</Link>";
  my $e=index($cd,$close,$from);
  my $blk=substr($cd,$from,$e-$from);
  my ($tag)=$blk=~/<Tag>([^<]*)<\/Tag>/;
  my ($txt)=$blk=~/<Text>([^<]*)<\/Text>/;
  push @ents,{kind=>$kind,tag=>$tag//"",text=>$txt//"",
              lead=>substr($open,0,index($open,'>')+1),blk=>$blk};
}
sub entry_close { $_[0]{kind} eq "N" ? "</Node>" : "</Link>" }

# ---------- 1+2: delete my appended link drawings + 4 orphans ----------
my %orphan = map { $_=>1 } (
  '8c69fe26-dab1-42dd-9019-bf8123a08789','40346d5f-8029-4655-b5bc-c74ec411af47',
  '5e1b390e-0362-4c3f-95ef-e0ca5579aa21','841fe248-4054-4d84-b7cc-6631d7128917');
my ($delMine,$delOrph)=(0,0);
my @keep;
for my $e (@ents) {
  if ($e->{kind} eq 'L' && $e->{tag} =~ /^cafe/) { $delMine++; next; }
  if ($e->{kind} eq 'L' && $orphan{$e->{tag}})   { $delOrph++; next; }
  push @keep,$e;
}
@ents=@keep;
print "deleted: myLinkDrawings=$delMine orphans=$delOrph\n";
die "expected 78 mine" unless $delMine==78;
die "expected 4 orphans" unless $delOrph==4;

# ---------- 3: retemplate my 8 adapter node drawings ----------
# original adapter drawing: find flat adapter node id from its Text
my $adTmplId;
{ my $t=index($x,'<Text>L Hip flx Ia Adapter</Text>');
  my $s=rindex($x,'<Node>',$t);
  ($adTmplId)=substr($x,$s,$t-$s)=~/<ID>([^<]+)<\/ID>/; }
die "adapter template node id\n" unless $adTmplId;
my ($adTmplEnt)=grep { $_->{kind} eq 'N' && $_->{tag} eq $adTmplId } @ents;
die "adapter drawing template not on page\n" unless $adTmplEnt;
my $retmpl=0;
for my $e (@ents) {
  next unless $e->{kind} eq 'N' && $e->{tag}=~/^cafe/ && $e->{text}=~/Adapter$/;
  my $new=$adTmplEnt->{blk};
  $new =~ s/<Text>[^<]*<\/Text>/<Text>$e->{text}<\/Text>/;
  $new =~ s/<Tag>[^<]*<\/Tag>/<Tag>$e->{tag}<\/Tag>/;
  $e->{blk}=$new; $e->{lead}=$adTmplEnt->{lead};
  $retmpl++;
}
print "retargeted adapter drawings: $retmpl\n";
die "expected 8" unless $retmpl==8;

# ---------- 4: add 4 contact-body node drawings ----------
my %bodyId;
for my $bn ("foot_L_contact","toe_L_contact","foot_R_contact","toe_R_contact") {
  my $i=index($x,"<Name>$bn</Name>");
  ($bodyId{$bn}) = substr($x,$i,300) =~ /<ID>([^<]+)<\/ID>/;
  die "body id $bn\n" unless $bodyId{$bn};
}
my $physTmplId='3c7d0a82-a0a4-4f8e-aed5-b40165b8c0bc';   # drawn muscle body
my ($physTmplEnt)=grep { $_->{kind} eq 'N' && $_->{tag} eq $physTmplId } @ents;
die "physical-body drawing template not found\n" unless $physTmplEnt;
for my $bn (sort keys %bodyId) {
  my $new=$physTmplEnt->{blk};
  $new =~ s/<Text>[^<]*<\/Text>/<Text>$bn<\/Text>/;
  $new =~ s/<Tag>[^<]*<\/Tag>/<Tag>$bodyId{$bn}<\/Tag>/;
  push @ents,{kind=>'N',tag=>$bodyId{$bn},text=>$bn,lead=>$physTmplEnt->{lead},blk=>$new};
}
print "added 4 contact-body drawings\n";

# ---------- 5: positions ----------
# occupancy from ALL node drawings
my @occ;
for my $e (@ents) {
  next unless $e->{kind} eq 'N';
  my ($l)=($e->{lead}.$e->{blk})=~/Left="([-\d.]+)"/; my ($t)=($e->{lead}.$e->{blk})=~/Top="([-\d.]+)"/;
  push @occ,[$l+$e->{blk}=~/Width="([\d.]+)"/?$1:32, $t] if defined $l && defined $t;
}
# recompute occupancy AFTER retarget (lead carries Left/Top in both formats)
@occ=[];
{
  @occ=();
  for my $e (@ents) {
    next unless $e->{kind} eq 'N';
    my $full=$e->{lead}.$e->{blk};
    my ($l)=$full=~/Left="([-\d.]+)"/; my ($t)=$full=~/Top="([-\d.]+)"/;
    next unless defined $l && defined $t;
    my ($w)=$full=~/Width="([\d.]+)"/; $w=32 unless $w;
    push @occ,[$l,$t,$l+$w,$t+40];
  }
}
sub free {
  my ($px,$py)=@_;
  for my $o (@occ) {   # $o=[x0,y0,x1,y1]
    return 0 if $px < $o->[2]+18 && $px+60 > $o->[0]-18 && $py < $o->[3]+14 && $py+40 > $o->[1]-14;
  }
  return 1;
}
sub claim { my ($px,$py)=@_; push @occ,[$px,$py,$px+60,$py+40]; }
sub place_near {   # find first free spot scanning right from (x,y)
  my ($px,$py)=@_;
  my $tries=0;
  while (!free($px,$py)) { $px+=46; $tries++; last if $tries>60; }
  claim($px,$py);
  return ($px,$py);
}

# helper: current position of a drawn node
sub posof { my ($tag)=@_;
  for my $e (@ents) { next unless $e->{kind} eq 'N' && $e->{tag} eq $tag;
    my $full=$e->{lead}.$e->{blk};
    my ($l)=$full=~/Left="([-\d.]+)"/; my ($t)=$full=~/Top="([-\d.]+)"/;
    return ($l,$t); }
  return (); }
sub setpos { my ($tag,$nx,$ny)=@_;
  for my $e (@ents) { next unless $e->{kind} eq 'N' && $e->{tag} =~ /^\Q$tag\E/;
    # Left/Top live in the opening tag (= lead)
    $e->{lead} =~ s/Left="[-\d.]+"/Left="$nx"/ or die "no Left on $tag\n";
    $e->{lead} =~ s/Top="[-\d.]+"/Top="$ny"/ or die "no Top on $tag\n";
    return 1; }
  die "setpos: drawing with tag $tag not found\n"; }

# reference nodes
my %REF;
for my $p (['LRGext','17fdb1f6'],['LRGextIN','6e188aac'],['LRGflx','7bd02f2b'],['LRGflxIN','136ad874'],
           [hipIaL=>'c92518bf'],[hipIaR=>'9cb6dfed'],[kneeIaL=>'7d681a6f'],[kneeIaR=>'4edb4e2b']) {
  my ($k,$pfx)=@$p;
  my $full=resolve($x,$pfx);
  my @p=posof($full);
  $REF{$k}=@p?[@p]:undef;
}
my @rgPlan = (
  ['cafe0001','R RG ext'],   ['cafe0002','R RG ext IN'],
  ['cafe0003','R RG flx'],   ['cafe0004','R RG flx IN']);
my @cmPlan = (['cafe0005','L c1'],['cafe0006','R c1'],['cafe0007','L V3'],['cafe0008','R V3']);
# R RG mirrors L RG block, shifted down 150
if ($REF{LRGext} && $REF{LRGextIN} && $REF{LRGflx} && $REF{LRGflxIN}) {
  my %mirror = (
    'cafe0001'=>[ @{$REF{LRGext}},   0,150 ],
    'cafe0002'=>[ @{$REF{LRGextIN}}, 0,150 ],
    'cafe0003'=>[ @{$REF{LRGflx}},   0,150 ],
    'cafe0004'=>[ @{$REF{LRGflxIN}}, 0,150 ]);
  for my $spec (@rgPlan) {
    my ($tag)=@$spec; my $m=$mirror{$tag};
    my ($nx,$ny)=place_near($m->[0]+$m->[2],$m->[1]+$m->[3]);
    setpos($tag,$nx,$ny); print "pos $spec->[1] -> ($nx,$ny)\n";
  }
} else {
  my ($nx,$ny);
  my $y0=760;
  for my $i (0..3) { ($nx,$ny)=place_near(180+($i%2)*90, $y0+int($i/2)*60); setpos($rgPlan[$i][0],$nx,$ny); }
}
# c1/V3 + II clusters: dedicated band right of everything (guaranteed free)
my ($maxX,$maxY)=(0,0);
for my $o (@occ) { $maxX=$o->[2] if $o->[2]>$maxX; $maxY=$o->[3] if $o->[3]>$maxY; }
{
  my $x0=$maxX+80; my $y=200;
  # c1/V3 2x2
  my @cm2=(['cafe0005','L c1'],['cafe0006','R c1'],['cafe0007','L V3'],['cafe0008','R V3']);
  for my $i (0..3) {
    my ($nx,$ny)=($x0+($i%2)*90, $y+int($i/2)*60);
    setpos($cm2[$i][0],$nx,$ny); claim($nx,$ny);
    print "pos $cm2[$i][1] -> ($nx,$ny)\n";
  }
  # II clusters: adapter + relay pairs, 4 rows
  my @ii2=(['cafe0111','cafe0101','L Hip flx II'],
           ['cafe0112','cafe0102','R Hip flx II'],
           ['cafe0113','cafe0103','L Knee flx II'],
           ['cafe0114','cafe0104','R Knee flx II']);
  for my $i (0..3) {
    my ($atag,$ntag,$nname)=@{$ii2[$i]};
    my ($ax,$ay)=($x0, $y+200+$i*80);
    my ($nx,$ny)=($x0+90, $ay);
    setpos($atag,$ax,$ay); claim($ax,$ay);
    setpos($ntag,$nx,$ny); claim($nx,$ny);
    print "pos $nname -> relay($nx,$ny), adapter($ax,$ay)\n";
  }
}
# contact groups at the bottom: body -> adapter -> neuron per side
{
  my $maxY=0; for my $o (@occ){ $maxY=$o->[3] if $o->[3]>$maxY; }
  my @rows=( ['foot_L_contact','cafe0211','L heel contact','cafe0201'],
             ['toe_L_contact', 'cafe0212','L toe contact', 'cafe0202'],
             ['foot_R_contact','cafe0213','R heel contact','cafe0203'],
             ['toe_R_contact', 'cafe0214','R toe contact', 'cafe0204'] );
  my $y=$maxY+70; my $x0=80;
  for my $i (0..3) {
    my ($bx,$by)=place_near($x0+$i*260,$y);
    setpos($bodyId{$rows[$i][0]},$bx,$by);
    my ($ax,$ay)=place_near($bx+90,$by);
    setpos($rows[$i][1],$ax,$ay);
    my ($nx,$ny)=place_near($ax+90,$by);
    setpos($rows[$i][3],$nx,$ny);
    print "pos $rows[$i][2] group -> body($bx,$by) adapter($ax,$ay) neuron($nx,$ny)\n";
  }
}

# ---------- 6: recompute indexes, normalize ALL drawings, re-add mine ----------
my %idx; my $i=0;
for my $e (@ents) { $idx{$e->{tag}}=$i++ if $e->{tag} && !exists $idx{$e->{tag}}; }
# normalize existing drawings (originals whose targets shifted past the
# deleted mid-file orphans, plus any other drift) to recomputed endpoints
my $norm=0;
for my $e (@ents) {
  next unless $e->{kind} eq 'L';
  next unless $e->{tag} && $func{$e->{tag}} && !($e->{tag}=~/^cafe/);
  my ($o,$d)=@{$func{$e->{tag}}};
  my ($oi,$di)=($idx{$o}//-1,$idx{$d}//-1);
  next if $oi<0 || $di<0;
  my $lead=$e->{lead};
  my ($oo)=$lead=~/Org="(\d+)"/; my ($dd)=$lead=~/Dst="(\d+)"/;
  next if defined $oo && $oo==$oi && defined $dd && $dd==$di;
  $lead =~ s/Org="\d+"/Org="$oi"/;
  $lead =~ s/Dst="\d+"/Dst="$di"/;
  $e->{lead}=$lead; $norm++;
}
print "normalized existing drawings: $norm\n";
# link template: an original synapse drawing WITHOUT Points (simple one)
my ($ltmpl);
for my $e (@ents) {
  next unless $e->{kind} eq 'L';
  next if $e->{blk}=~/<Point /;
  next unless $e->{tag} && exists $func{$e->{tag}};
  $ltmpl=$e->{blk}; last;
}
$ltmpl or die "no simple link drawing template\n";
my @myLinks = grep { /^cafe/ } sort keys %func;
die "expected 78 cafe functional links, got ".scalar(@myLinks)."\n" unless @myLinks==78;
my $added=0; my $missing=0;
for my $lid (@myLinks) {
  my ($o,$d)=@{$func{$lid}};
  my ($oi,$di)=($idx{$o}//-1,$idx{$d}//-1);
  if ($oi<0 || $di<0) { print "WARN: endpoint not drawn for $lid ($o/$d)\n"; $missing++; next; }
  my $n=$ltmpl;
  $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$lid<\/Tag>/;
  $n =~ s/Org="\d+"/Org="$oi"/;
  $n =~ s/Dst="\d+"/Dst="$di"/;
  $n =~ s/\s*<Point [^>]*\/>//g;
  push @ents,{kind=>'L',tag=>$lid,text=>'',lead=>'<Link Org="'.$oi.'" Dst="'.$di.'">',blk=>$n};
  $added++;
}
print "re-added link drawings: $added (missing endpoints: $missing)\n";

# ---------- 7: rebuild CDATA ----------
my $out='<Root>
<Diagram>
';
# preserve everything up to <AddFlow ...> from original cd header
my ($hdr) = $cd =~ /^(\s*<Root>.*?<AddFlow Nodes="\d+" Links="\d+">\s*)/s
  or die "no header\n";
my $nn = scalar(grep{$_->{kind} eq 'N'}@ents);
my $ln2 = scalar(grep{$_->{kind} eq 'L'}@ents);
$hdr =~ s/Nodes="\d+"/Nodes="$nn"/;
$hdr =~ s/Links="\d+"/Links="$ln2"/;
my $body='';
for my $e (@ents) {
  $body .= '  '.$e->{lead}."\n".$e->{blk}."\n".entry_close($e)."\n";
}
my $newcd = $hdr.$body.'  </AddFlow>
</Diagram>
</Root>';
substr($x,$CS,$CE-$CS) = $newcd;
spew($PROJ,$x);
my $nn=scalar(grep{$_->{kind} eq 'N'}@ents); my $ln=scalar(grep{$_->{kind} eq 'L'}@ents);
print "page: Nodes=$nn Links=$ln. written.\n";
