#!/usr/bin/perl
# build_contact.pl — Walker_2_Layer_CPG_BilateralRG: heel/toe contact drive.
# Per side: "L heel contact" + "L toe contact" nonspiking neurons, fed by
# PhysicalToNodeAdapters (SourceDataType=ContactCount, TargetDataType=
# ExternalCurrent, Gain C=20) from foot_L_contact / toe_L_contact bodies,
# each sending EXCITATORY synapses ("Afferent HC Excite" type, weak) to the
# ipsilateral extensor layers: RG ext, Hip PF ext, Knee PF ext, Hip MN ext,
# Knee MN ext. 4 neurons, 4 adapters (+8 adapter links), 20 synapses.
# Usage: perl build_contact.pl <dir>
use strict; use warnings;
my $dir = shift @ARGV or die "usage: build_contact.pl <dir>\n";
my $PROJ = "$dir/Walker_2_Layer_CPG_BilateralRG.aproj";
my $ASIM = "$dir/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim";
sub slurp { open my $f,'<',$_[0] or die "$_[0]: $!"; local $/; my $s=<$f>; close $f; $s }
sub spew  { open my $o,'>',$_[0] or die "$_[0]: $!"; print {$o} $_[1]; close $o }
sub g { my $n=shift; sprintf("cafe%04d-0000-4000-8000-%012d",$n,900000000000+$n) }
my $childseq = 1300;
sub childguid { g($childseq++) }
sub reroll_ids { my ($blk,$main)=@_; my $first=1;
  $blk =~ s{<ID>([^<]+)</ID>}{ $first ? do{$first=0;"<ID>$main</ID>"} : "<ID>".childguid()."</ID>" }ge; $blk }
sub resolve { my ($x,$p)=@_; my $i=index($x,"<ID>$p"); die "prefix $p not found\n" if $i<0;
  substr($x,$i+4,index($x,'</ID>',$i+4)-$i-4) }

# body ids resolved by NAME (contact bodies)
my %P = (
  rgExtL=>'17fdb1f6', rgExtR=>'cafe0001',
  pfHipEL=>'57b344bf', pfHipER=>'aa4d0055', pfKneeEL=>'f336ecdc', pfKneeER=>'eb127c29',
  mnHipEL=>'64083a38', mnHipER=>'56e26621', mnKneeEL=>'c0986b70', mnKneeER=>'cabac9d1',
  affType=>'cafe0043',
);
my @NEUR = ( ['heelL','L heel contact','foot_L_contact'], ['toeL','L toe contact','toe_L_contact'],
             ['heelR','R heel contact','foot_R_contact'], ['toeR','R toe contact','toe_R_contact'] );
my %N = ( heelL=>g(201), toeL=>g(202), heelR=>g(203), toeR=>g(204) );
my @AD = (g(211),g(212),g(213),g(214));
my @ADL = ((map { g(221+$_) } 0..7));
my @LK = ((map { g(231+$_) } 0..19));
# targets per side: rgExt, pfHipE, pfKneeE, mnHipE, mnKneeE
my %T = ( L=>[qw(rgExtL pfHipEL pfKneeEL mnHipEL mnKneeEL)], R=>[qw(rgExtR pfHipER pfKneeER mnHipER mnKneeER)] );

# ---- ASIM ----
{
  my $x = slurp($ASIM);
  if (index($x,"<Name>L heel contact</Name>")>=0) { print "asim: contact layer present, skip\n"; }
  else {
  my %R = ( (map { $_=>resolve($x,$P{$_}) } keys %P), %N );
  # body ids by name
  my %body;
  for my $bn ("foot_L_contact","toe_L_contact","foot_R_contact","toe_R_contact") {
    my $i=index($x,"<Name>$bn</Name>"); die "body $bn\n" if $i<0;
    my $s=rindex($x,"<RigidBody>",$i);
    ($body{$bn}) = substr($x,$i,300) =~ /<ID>([^<]+)<\/ID>/;
    die "body id $bn\n" unless $body{$bn};
  }
  # contact neurons: clone an existing plain afferent relay ("L Hip flx II" or the Ia relay)
  { my $tid=resolve($x,'c92518bf');
    my $i=index($x,"<ID>$tid</ID>"); die "tmpl\n" if $i<0;
    my $s=rindex($x,"<Neuron>",$i); my $e=index($x,"</Neuron>",$i)+9;
    my $blk=substr($x,$s,$e-$s); my $at=$e;
    for my $spec (@NEUR) {
      my ($key,$name,$bodyname)=@$spec;
      my $new=reroll_ids($blk,$N{$key});
      $new =~ s/<Name>[^<]*<\/Name>/<Name>$name<\/Name>/ or die;
      substr($x,$at,0)="\n".$new; $at+=length($new)+1; } }
  # adapters: clone the asim II adapter (PhysicalToNode) and repoint src=contact body, DataType=ContactCount, Gain C
  { my $ti=index($x,"<Name>L Hip flx II Adapter</Name>"); die "II adpt tmpl\n" if $ti<0;
    my $s=rindex($x,"<Adapter>",$ti); my $e=index($x,"</Adapter>",$ti)+10;
    my $blk=substr($x,$s,$e-$s); my $at=$e;
    for my $k (0..$#NEUR) {
      my ($key,$name,$bodyname)=@{$NEUR[$k]};
      my $new=reroll_ids($blk,$AD[$k]);
      my ($osrc)=$new=~/<SourceID>([^<]+)<\/SourceID>/;
      my ($odst)=$new=~/<TargetID>([^<]+)<\/TargetID>/;
      $new =~ s/<Name>[^<]*<\/Name>/<Name>$name Adapter<\/Name>/ or die;
      $new =~ s/<SourceID>\Q$osrc\E<\/SourceID>/<SourceID>$body{$bodyname}<\/SourceID>/ or die;
      $new =~ s/<TargetID>\Q$odst\E<\/TargetID>/<TargetID>$N{$key}<\/TargetID>/ or die;
      $new =~ s/<SourceDataType>[^<]*<\/SourceDataType>/<SourceDataType>ContactCount<\/SourceDataType>/ or die;
      $new =~ s/(<C>)[^<]*(<\/C>)/${1}20$2/ or print "WARN: no <C> gain in adapter clone (gain untouched)\n";
      substr($x,$at,0)="\n".$new; $at+=length($new)+1; } }
  # 20 synapses
  { my ($to,$td)=($R{rgExtL},$R{mnHipEL});
    # template: any existing Afferent HC Excite connexion: find one by type id
    my ($found); pos($x)=0;
    while ($x =~ /<Connexion>(.*?)<\/Connexion>/gs) {
      my $b=$1; my ($sy)=$b=~/<SynapseTypeID>([^<]+)/;
      next unless $sy && $sy eq $R{affType};
      $found=$b; ($to,$td)=($b=~/<SourceID>([^<]+)/,$b=~/<TargetID>([^<]+)/); last; }
    die "no Afferent HC Excite connexion template\n" unless $found;
    my $blk="<Connexion>$found</Connexion>";
    my $at=index($x,$blk)+length($blk); my $ins=''; my $li=0;
    for my $side (qw(L R)) {
      for my $key ($side eq 'L' ? qw(heelL toeL) : qw(heelR toeR)) {
        for my $tp (@{$T{$side}}) {
          my $new=reroll_ids($blk,$LK[$li++]);
          $new =~ s/<SourceID>\Q$to\E<\/SourceID>/<SourceID>$N{$key}<\/SourceID>/ or die;
          $new =~ s/<TargetID>\Q$td\E<\/TargetID>/<TargetID>$R{$tp}<\/TargetID>/ or die;
          $ins .= $new."\n"; } } }
    substr($x,$at,0)="\n".$ins; }
  spew($ASIM,$x); print "asim written\n";
  }
}

# ---- APROJ ----
{
  my $x = slurp($PROJ);
  die "aproj already has contact layer\n" if $x=~/<Text>L heel contact<\/Text>/;
  my %R = ( (map { $_=>resolve($x,$P{$_}) } keys %P), %N );
  my %body;
  for my $bn ("foot_L_contact","toe_L_contact","foot_R_contact","toe_R_contact") {
    my $i=index($x,"<Name>$bn</Name>"); die "body $bn\n" if $i<0;
    ($body{$bn}) = substr($x,$i,300) =~ /<ID>([^<]+)<\/ID>/;
    die "body id $bn\n" unless $body{$bn};
  }
  sub node_block { my ($x,$id)=@_; my $i=index($x,"<ID>$id</ID>"); return () if $i<0;
    my $s=rindex($x,'<Node>',$i); my $e=index($x,'</Node>',$i)+7; ($s<0||$e<7)?():($s,$e-$s) }
  sub add_to_list { my ($xr,$nid,$lid,$list)=@_;
    my ($s,$l)=node_block($$xr,$nid); die "node $nid list add\n" unless $s;
    my $seg=substr($$xr,$s,$l); my $ns=$seg;
    my $rep="<$list>\n<ID>$lid</ID>\n</$list>";
    if ($ns=~s{<$list/>\s*}{$rep}s){} elsif ($ns=~s{<$list>\s*</$list>}{$rep}s){}
    elsif ($ns=~s{</$list>}{<ID>$lid</ID>\n</$list>}s){} else {die "no $list on $nid\n"}
    substr($$xr,$s,$l)=$ns }
  sub find_adpt_link { my ($x,$o,$d)=@_; pos($$x)=0;
    while ($$x =~ /<Link>(.*?)<\/Link>/gs) { my $b=$1;
      next unless $b =~ /Behavior\.Links\.Adapter</;
      my ($bo)=$b=~/<OriginID>([^<]+)/; my ($bd)=$b=~/<DestinationID>([^<]+)/;
      next unless defined $bo && defined $bd;
      return "<Link>$b</Link>" if $bo eq $o && $bd eq $d }
    undef }

  # contact neurons: clone the "L Hip flx II" node
  { my ($s,$l)=node_block($x,resolve($x,'cafe0101')); die "II node tmpl\n" unless $s;
    my $blk=substr($x,$s,$l);
    for my $spec (@NEUR) {
      my ($key,$name,$bn)=@$spec;
      my $new=reroll_ids($blk,$N{$key});
      $new =~ s/<Text>[^<]*<\/Text>/<Text>$name<\/Text>/ or die;
      $new =~ s{<InLinks>.*?</InLinks>}{<InLinks>\n</InLinks>}gs;
      $new =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks>\n</OutLinks>}gs;
      substr($x,$s+$l,0)="\n".$new; } }
  # adapter nodes: clone "L Hip flx II Adapter", repoint OriginID=contact body, DataTypeID=ContactCount
  { my $at=index($x,"<Text>L Hip flx II Adapter</Text>"); die "adpt tmpl\n" if $at<0;
    my $s=rindex($x,"<Node>",$at); my $e=index($x,"</Node>",$at)+7;
    my $blk=substr($x,$s,$e-$s);
    my ($iaAdId)=$blk=~/<ID>([^<]+)<\/ID>/;
    for my $k (0..$#NEUR) {
      my ($key,$name,$bn)=@{$NEUR[$k]};
      my $new=reroll_ids($blk,$AD[$k]);
      my ($odst)=$new=~/<DestinationID>([^<]+)<\/DestinationID>/;
      $new =~ s/<Text>[^<]*<\/Text>/<Text>$name Adapter<\/Text>/ or die;
      $new =~ s/<OriginID>[^<]+<\/OriginID>/<OriginID>$body{$bn}<\/OriginID>/ or die;
      $new =~ s/<DestinationID>\Q$odst\E<\/DestinationID>/<DestinationID>$N{$key}<\/DestinationID>/ or die;
      $new =~ s/<DataTypeID>[^<]*<\/DataTypeID>/<DataTypeID>ContactCount<\/DataTypeID>/ or die;
      $new =~ s/(<C Value=")[^"]*(" Scale="nano"[^>]*>)/${1}20$2/ or print "WARN aproj gain C untouched\n";
      $new =~ s{<InLinks>.*?</InLinks>}{<InLinks>\n</InLinks>}gs;
      $new =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks>\n</OutLinks>}gs;
      substr($x,$e,0)="\n".$new; }
    # adapter links: receptor->adpt tmpl = (IIsrc->IIadpt); adpt->node tmpl = (IIadpt->IInode)
    my ($iisrc)= $blk =~ /<OriginID>([^<]+)<\/OriginID>/;
    my ($iinodedst) = $blk =~ /<DestinationID>([^<]+)<\/DestinationID>/;
    my $al1=find_adpt_link(\$x,$iisrc,$iaAdId);
    my $al2=find_adpt_link(\$x,$iaAdId,$iinodedst);
    die "adapter link templates missing\n" unless $al1 && $al2;
    for my $k (0..$#NEUR) {
      my ($key,$name,$bn)=@{$NEUR[$k]};
      my $a1=reroll_ids($al1,$ADL[$k*2]);
      $a1 =~ s/<OriginID>[^<]+<\/OriginID>/<OriginID>$body{$bn}<\/OriginID>/ or die;
      $a1 =~ s/<DestinationID>[^<]+<\/DestinationID>/<DestinationID>$AD[$k]<\/DestinationID>/ or die;
      my $a2=reroll_ids($al2,$ADL[$k*2+1]);
      $a2 =~ s/<OriginID>[^<]+<\/OriginID>/<OriginID>$AD[$k]<\/OriginID>/ or die;
      $a2 =~ s/<DestinationID>[^<]+<\/DestinationID>/<DestinationID>$N{$key}<\/DestinationID>/ or die;
      substr($x,index($x,$al2)+length($al2),0)="\n".$a1."\n".$a2;
      # body sources have no neural Node block -> no OutLinks registration
      add_to_list(\$x,$AD[$k],$ADL[$k*2],'InLinks');
      add_to_list(\$x,$AD[$k],$ADL[$k*2+1],'OutLinks');
      add_to_list(\$x,$N{$key},$ADL[$k*2+1],'InLinks'); } }
  # 20 synapse links: template = first afferent link (type cafe0043) — use the iaHipL->pfHipFL link
  { my ($t); pos($x)=0;
    while ($x =~ /<Link>(.*?)<\/Link>/gs) { my $b=$1;
      next unless $b =~ /Behavior\.Synapse</;
      my ($sy)=$b=~/<SynapticTypeID>([^<]+)/;
      next unless $sy && $sy eq $R{affType};
      $t="<Link>$b</Link>"; last; }
    die "aff synapse tmpl missing\n" unless $t;
    my ($to,$td) = ($t =~ /<OriginID>([^<]+)/, $t =~ /<DestinationID>([^<]+)/);
    my $at=index($x,$t)+length($t); my $ins=''; my $li=0;
    for my $side (qw(L R)) {
      for my $key ($side eq 'L' ? qw(heelL toeL) : qw(heelR toeR)) {
        for my $tp (@{$T{$side}}) {
          my $new=reroll_ids($t,$LK[$li++]);
          $new =~ s/<OriginID>\Q$to\E<\/OriginID>/<OriginID>$N{$key}<\/OriginID>/ or die;
          $new =~ s/<DestinationID>\Q$td\E<\/DestinationID>/<DestinationID>$R{$tp}<\/DestinationID>/ or die;
          $ins .= $new."\n"; } } }
    substr($x,$at,0)="\n".$ins;
    $li=0;
    for my $side (qw(L R)) {
      for my $key ($side eq 'L' ? qw(heelL toeL) : qw(heelR toeR)) {
        for my $tp (@{$T{$side}}) {
          add_to_list(\$x,$N{$key},$LK[$li],'OutLinks');
          add_to_list(\$x,$R{$tp},$LK[$li],'InLinks'); $li++; } } }
    # drawings
    my ($cd)=$x =~ /<DiagramXml><!\[CDATA\[(.*?)\]\]><\/DiagramXml>/s;
    my ($tmplnode)=$cd =~ m{(<Node Left="[^"]+" Top="[^"]+"[^>]*>(?:(?!</Node>).)*?<Tag>\Q$R{rgExtL}\E</Tag>(?:(?!</Node>).)*?</Node>)}s;
    die "drawing tmpl\n" unless $tmplnode;
    my $ins2='';
    for my $k (0..$#NEUR) {
      my $n=$tmplnode; my $top=1000+$k*40;
      $n =~ s/Left="[-\d.]+"/Left="420"/; $n =~ s/Top="[-\d.]+"/Top="$top"/;
      $n =~ s/<Text>[^<]*<\/Text>/<Text>$NEUR[$k][1]<\/Text>/;
      $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$N{$NEUR[$k][0]}<\/Tag>/;
      $ins2 .= "  ".$n."\n";
      my $a=$tmplnode;
      $a =~ s/Left="[-\d.]+"/Left="490"/; $a =~ s/Top="[-\d.]+"/Top="$top"/;
      $a =~ s/<Text>[^<]*<\/Text>/<Text>$NEUR[$k][1] Adapter<\/Text>/;
      $a =~ s/<Tag>[^<]*<\/Tag>/<Tag>$AD[$k]<\/Tag>/;
      $ins2 .= "  ".$a."\n"; }
    my ($tib)=$t=~/<ID>([^<]+)/;
    my ($tl)=$cd =~ m{(<Link Org="\d+" Dst="\d+">(?:(?!</Link>).)*?<Tag>\Q$tib\E</Tag>(?:(?!</Link>).)*?</Link>)}s;
    for my $i (0..19) { my $n=$tl; $n=~s/<Tag>[^<]*<\/Tag>/<Tag>$LK[$i]<\/Tag>/; $ins2 .= "  $n\n"; }
    for my $i (0..7) { my $n=$tl; $n=~s/<Tag>[^<]*<\/Tag>/<Tag>$ADL[$i]<\/Tag>/; $ins2 .= "  $n\n"; }
    $cd =~ s{(\s*</AddFlow>)}{\n$ins2$1}s;
    my $nn=()=$cd=~/<Node Left=/g; my $ln=()=$cd=~/<Link Org=/g;
    $cd =~ s{<AddFlow Nodes="\d+" Links="\d+"}{<AddFlow Nodes="$nn" Links="$ln"};
    $x =~ s{(<DiagramXml><!\[CDATA\[).*?(\]\]></DiagramXml>)}{$1$cd$2}s;
    print "aproj page: Nodes=$nn Links=$ln\n"; }
  spew($PROJ,$x); print "aproj written\n";
}
print "BUILD CONTACT DONE\n";
