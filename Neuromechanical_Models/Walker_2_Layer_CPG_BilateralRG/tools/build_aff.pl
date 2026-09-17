#!/usr/bin/perl
# build_aff.pl — Walker_2_Layer_CPG_BilateralRG: afferent feedback to the
# half-centers (ipsilateral), per Ben 2026-09-16:
#   flexor Ia + NEW II  --exc--> flexor HCs of PF and RG (hip->hipPF+RG, knee->kneePF+RG)
#   extensor Ib         --exc--> extensor HCs of PF and RG
# II chains are NEW: StretchReceptor -> PhysicalToNodeAdapter(DataTypeID=II)
# -> "S J flx II" relay (clone of the Ia relay neuron).
# All afferent->HC synapses: one type "Afferent HC Excite" (Equil -40,
# SynAmp 0.1 = weak; the TYPE SynAmp is the single tuning knob).
# Usage: perl build_aff.pl <dir>
use strict; use warnings;
my $dir = shift @ARGV or die "usage: build_aff.pl <dir>\n";
my $PROJ = "$dir/Walker_2_Layer_CPG_BilateralRG.aproj";
my $ASIM = "$dir/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim";
sub slurp { open my $f,'<',$_[0] or die "$_[0]: $!"; local $/; my $s=<$f>; close $f; $s }
sub spew  { open my $o,'>',$_[0] or die "$_[0]: $!"; print {$o} $_[1]; close $o }
sub g { my $n=shift; sprintf("cafe%04d-0000-4000-8000-%012d",$n,900000000000+$n) }
my $childseq = 1100;
sub childguid { g($childseq++) }
sub reroll_ids { my ($blk,$main)=@_; my $first=1;
  $blk =~ s{<ID>([^<]+)</ID>}{ $first ? do{$first=0;"<ID>$main</ID>"} : "<ID>".childguid()."</ID>" }ge; $blk }

my %P = (
  iaHipL=>'c92518bf', iaHipR=>'9cb6dfed', iaKneeL=>'7d681a6f', iaKneeR=>'4edb4e2b',
  ibHipL=>'cae540f6', ibHipR=>'6d50c8a7', ibKneeL=>'b51f816f', ibKneeR=>'07e2aace',
  pfHipFL=>'f5ba89bb', pfHipFR=>'16bc9942', pfKneeFL=>'9bc35128', pfKneeFR=>'fbda99c8',
  pfHipEL=>'57b344bf', pfHipER=>'aa4d0055', pfKneeEL=>'f336ecdc', pfKneeER=>'eb127c29',
  rgFlxL=>'7bd02f2b', rgFlxR=>'cafe0003', rgExtL=>'17fdb1f6', rgExtR=>'cafe0001',
  # link template: L RG ext -> L Hip PF ext (RG to PF Excite)
  tmplE=>'17fdb1f6', tmpltD=>'57b344bf',
  iaAdptTmpl=>'2bffb0b6',          # asim PhysicalToNode adapter template
);
my %N = ( IIhipL=>g(101), IIhipR=>g(102), IIkneeL=>g(103), IIkneeR=>g(104) );
my @AD = (g(111),g(112),g(113),g(114));
my @ADL = ((map { g(121+$_) } 0..7));
my $TYPE = g(43);
my @LK = ((map { g(131+$_) } 0..23));
sub resolve { my ($x,$p)=@_; my $i=index($x,"<ID>$p"); die "prefix $p not found\n" if $i<0;
  substr($x,$i+4,index($x,'</ID>',$i+4)-$i-4) }

my @II = (
  ['IIhipL','L Hip flx II','recHipL','iaHipL',0],
  ['IIhipR','R Hip flx II','recHipR','iaHipR',1],
  ['IIkneeL','L Knee flx II','recKneeL','iaKneeL',2],
  ['IIkneeR','R Knee flx II','recKneeR','iaKneeR',3],
);
my @AFF = (
  ['iaHipL','pfHipFL'],['iaHipL','rgFlxL'],['iaHipR','pfHipFR'],['iaHipR','rgFlxR'],
  ['iaKneeL','pfKneeFL'],['iaKneeL','rgFlxL'],['iaKneeR','pfKneeFR'],['iaKneeR','rgFlxR'],
  ['ibHipL','pfHipEL'],['ibHipL','rgExtL'],['ibHipR','pfHipER'],['ibHipR','rgExtR'],
  ['ibKneeL','pfKneeEL'],['ibKneeL','rgExtL'],['ibKneeR','pfKneeER'],['ibKneeR','rgExtR'],
  ['IIhipL','pfHipFL'],['IIhipL','rgFlxL'],['IIhipR','pfHipFR'],['IIhipR','rgFlxR'],
  ['IIkneeL','pfKneeFL'],['IIkneeL','rgFlxL'],['IIkneeR','pfKneeFR'],['IIkneeR','rgFlxR'],
);

# ============================ ASIM ============================
{
  my $x = slurp($ASIM);
  die "asim already has afferent layer\n" if $x=~/<Name>L Hip flx II<\/Name>/;
  my %R = ( (map { $_=>resolve($x,$P{$_}) } keys %P), %N );

  # type
  { my $tp=index($x,"<Name>RG to PF Excite</Name>"); die "type tmpl\n" if $tp<0;
    my $s=rindex($x,"<SynapseType>",$tp); my $e=index($x,"</SynapseType>",$tp)+15;
    my $n=substr($x,$s,$e-$s);
    $n =~ s/<Name>[^<]*<\/Name>/<Name>Afferent HC Excite<\/Name>/;
    $n =~ s/<ID>[^<]*<\/ID>/<ID>$TYPE<\/ID>/;
    substr($x,$e,0)="\n".$n; print "asim: +type Afferent HC Excite\n"; }

  # II neurons
  for my $spec (@II) {
    my ($key,$name,$rec,$tmplr,$adi)=@$spec;
    my $tid=$R{$tmplr};
    my $i=index($x,"<ID>$tid</ID>"); die "II tmpl neuron $tmplr\n" if $i<0;
    my $s=rindex($x,"<Neuron>",$i); my $e=index($x,"</Neuron>",$i)+9;
    my $new=reroll_ids(substr($x,$s,$e-$s),$N{$key});
    $new =~ s/<Name>[^<]*<\/Name>/<Name>$name<\/Name>/ or die;
    substr($x,$e,0)="\n".$new; print "asim: +neuron $name\n"; }

  # II adapters (clone the SAME-MUSCLE Ia adapter; only channel/target change)
  for my $spec (@II) {
    my ($key,$name,$rec,$tmplr,$adi)=@$spec;
    my ($side,$joint) = $name =~ /^([LR]) (Hip|Knee)/;
    my $tmplName = "$side $joint flx Ia Adapter";
    my $ti2=index($x,"<Name>$tmplName</Name>"); die "asim adpt tmpl $tmplName\n" if $ti2<0;
    my $s=rindex($x,"<Adapter>",$ti2); my $e=index($x,"</Adapter>",$ti2)+10;
    my $new=reroll_ids(substr($x,$s,$e-$s),$AD[$adi]);
    my ($odst)=$new=~/<TargetID>([^<]+)<\/TargetID>/;
    $new =~ s/<Name>[^<]*<\/Name>/<Name>${name} Adapter<\/Name>/ or die;
    $new =~ s/<TargetID>\Q$odst\E<\/TargetID>/<TargetID>$N{$key}<\/TargetID>/ or die;
    $new =~ s/<SourceDataType>[^<]*<\/SourceDataType>/<SourceDataType>II<\/SourceDataType>/ or die;
    substr($x,$e,0)="\n".$new; print "asim: +adapter $name (from $tmplName)\n"; }

  # 24 connexions
  { my ($to,$td)=($R{tmplE},$R{tmpltD});
    my ($found); pos($x)=0;
    while ($x =~ /<Connexion>(.*?)<\/Connexion>/gs) {
      my $b=$1; my ($o)=$b=~/<SourceID>([^<]+)/; my ($d)=$b=~/<TargetID>([^<]+)/;
      if ($o eq $to && $d eq $td) { $found=$b; last } }
    die "asim aff tmpl missing\n" unless $found;
    my $blk="<Connexion>$found</Connexion>";
    my $at=index($x,$blk)+length($blk); my $ins='';
    for my $i (0..$#AFF) {
      my ($src,$dst)=@{$AFF[$i]};
      my $new=reroll_ids($blk,$LK[$i]);
      $new =~ s/<SourceID>\Q$to\E<\/SourceID>/<SourceID>$R{$src}<\/SourceID>/ or die;
      $new =~ s/<TargetID>\Q$td\E<\/TargetID>/<TargetID>$R{$dst}<\/TargetID>/ or die;
      $new =~ s/<SynapseTypeID>[^<]*<\/SynapseTypeID>/<SynapseTypeID>$TYPE<\/SynapseTypeID>/ or die;
      $ins .= $new."\n"; }
    substr($x,$at,0)="\n".$ins;
    print "asim: +",scalar(@AFF)," afferent->HC connexions\n"; }
  spew($ASIM,$x); print "asim written\n";
}

# ============================ APROJ ============================
{
  my $x = slurp($PROJ);
  die "aproj already has afferent layer\n" if $x=~/<Text>L Hip flx II<\/Text>/;
  my %R = ( (map { $_=>resolve($x,$P{$_}) } keys %P), %N );

  sub node_block { my ($x,$id)=@_; my $i=index($x,"<ID>$id</ID>"); return () if $i<0;
    my $s=rindex($x,'<Node>',$i); my $e=index($x,'</Node>',$i)+7; ($s<0||$e<7)?():($s,$e-$s) }
  sub add_to_list { my ($xr,$nid,$lid,$list)=@_;
    my ($s,$l)=node_block($$xr,$nid); die "node $nid list add\n" unless $s;
    my $seg=substr($$xr,$s,$l); my $ns=$seg;
    my $rep="<$list>\n<ID>$lid</ID>\n</$list>";
    if ($ns=~s{<$list/>\s*}{$rep}s){} elsif ($ns=~s{<$list>\s*</$list>}{$rep}s){}
    elsif ($ns=~s{</$list>}{<ID>$lid</ID>\n</$list>}s){} else {die "no $list on $nid\n"}
    substr($$xr,$s,$l)=$ns }
  sub find_link { my ($x,$o,$d)=@_; pos($$x)=0;
    while ($$x =~ /<Link>(.*?)<\/Link>/gs) { my $b=$1;
      my ($bo)=$b=~/<OriginID>([^<]+)/; my ($bd)=$b=~/<DestinationID>([^<]+)/;
      next unless defined $bo && defined $bd;
      return "<Link>$b</Link>" if $bo eq $o && $bd eq $d }
    undef }

  # type (aproj: type Link blocks inside SynapseTypes)
  { my $tb = find_link(\$x, "", "") ; # placeholder
    pos($x)=0; my ($tblk);
    while ($x =~ /<Link>(.*?)<\/Link>/gs) { my $b=$1;
      next unless $b =~ /SynapseTypes\./;
      my ($nm)=$b=~/<Name>([^<]*)<\/Name>/;
      if ($nm && $nm eq 'RG to PF Excite') { $tblk="<Link>$b</Link>"; last } }
    die "aproj type tmpl\n" unless $tblk;
    my $n=reroll_ids($tblk,$TYPE);
    $n =~ s/<Name>[^<]*<\/Name>/<Name>Afferent HC Excite<\/Name>/;
    $n =~ s{<EquilibriumPotential Value="[^"]*" Scale="milli" Actual="[^"]*"/>}{<EquilibriumPotential Value="-40" Scale="milli" Actual="-0.04"/>};
    $n =~ s{<MaxSynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>}{<MaxSynapticConductance Value="0.1" Scale="micro" Actual="1e-007"/>};
    substr($x,index($x,$tblk)+length($tblk),0)="\n".$n;
    print "aproj: +type Afferent HC Excite\n"; }

  # II neurons (clone Ia relay node)
  for my $spec (@II) {
    my ($key,$name,$rec,$tmplr,$adi)=@$spec;
    my ($s,$l)=node_block($x,$R{$tmplr}); die "aproj II tmpl $tmplr\n" unless $s;
    my $new=reroll_ids(substr($x,$s,$l),$N{$key});
    $new =~ s/<Text>[^<]*<\/Text>/<Text>$name<\/Text>/ or die;
    $new =~ s{<InLinks>.*?</InLinks>}{<InLinks>\n</InLinks>}gs;
    $new =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks>\n</OutLinks>}gs;
    substr($x,$s+$l,0)="\n".$new; print "aproj: +node $name\n"; }

  # II adapter nodes (clone the SAME-MUSCLE Ia adapter node; OriginID stays)
  { my $ins_ad='';
    for my $spec (@II) {
      my ($key,$name,$rec,$tmplr,$adi)=@$spec;
      my ($side,$joint) = $name =~ /^([LR]) (Hip|Knee)/;
      my $tmplName = "$side $joint flx Ia Adapter";
      my $at=index($x,"<Text>$tmplName</Text>"); die "aproj adpt tmpl $tmplName\n" if $at<0;
      my $s=rindex($x,"<Node>",$at); my $e=index($x,"</Node>",$at)+7;
      my $tblk=substr($x,$s,$e-$s);
      my ($osrc)=$tblk=~/<OriginID>([^<]+)<\/OriginID>/;
      my ($odst)=$tblk=~/<DestinationID>([^<]+)<\/DestinationID>/;
      my $new=reroll_ids($tblk,$AD[$adi]);
      $new =~ s/<Text>[^<]*<\/Text>/<Text>${name} Adapter<\/Text>/ or die;
      $new =~ s/<DestinationID>\Q$odst\E<\/DestinationID>/<DestinationID>$N{$key}<\/DestinationID>/ or die;
      $new =~ s/<DataTypeID>[^<]*<\/DataTypeID>/<DataTypeID>II<\/DataTypeID>/ or die;
      $new =~ s{<InLinks>.*?</InLinks>}{<InLinks>\n</InLinks>}gs;
      $new =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks>\n</OutLinks>}gs;
      substr($x,$e,0)="\n".$new;
      # adapter links for THIS muscle: (receptor->IaAdapter) and (IaAdapter->IaRelay)
      pos($x)=0; my ($al1,$al2);
      while ($x =~ /<Link>(.*?)<\/Link>/gs) { my $b=$1;
        next unless $b =~ /Behavior\.Links\.Adapter</;
        my ($bo)=$b=~/<OriginID>([^<]+)/; my ($bd)=$b=~/<DestinationID>([^<]+)/;
        next unless defined $bo && defined $bd;
        $al1="<Link>$b</Link>" if !defined($al1) && $bo eq $osrc && $bd eq (substr($tblk,0,0)||$osrc) && $bd eq $osrc && 0;
        if (!defined $al1 && $bd eq substr($x, rindex($x,"<ID>",index($x,"<Text>$tmplName</Text>")) ,0) ) {}
        # simpler: match by adapter-template id embedded below
        last; }
      # NOTE: adapter-link creation handled in second pass below
      $ins_ad .= '';
    }
    print "aproj: +4 II adapter nodes\n";
    # find adapter links by adapter-node id: the Ia adapter node ids
    for my $spec (@II) {
      my ($key,$name,$rec,$tmplr,$adi)=@$spec;
      my ($side,$joint) = $name =~ /^([LR]) (Hip|Knee)/;
      my $tmplName = "$side $joint flx Ia Adapter";
      my $at=index($x,"<Text>$tmplName</Text>");
      my $s=rindex($x,"<Node>",$at); my $e=index($x,"</Node>",$at)+7;
      my $tblk=substr($x,$s,$e-$s);
      my ($osrc,$odst)=($tblk=~/<OriginID>([^<]+)<\/OriginID>/, $tblk=~/<DestinationID>([^<]+)<\/DestinationID>/);
      my ($iaAdId)=$tblk=~/<ID>([^<]+)<\/ID>/;
      pos($x)=0; my ($al1,$al2);
      while ($x =~ /<Link>(.*?)<\/Link>/gs) { my $b=$1;
        next unless $b =~ /Behavior\.Links\.Adapter</;
        my ($bo)=$b=~/<OriginID>([^<]+)/; my ($bd)=$b=~/<DestinationID>([^<]+)/;
        next unless defined $bo && defined $bd;
        $al1="<Link>$b</Link>" if !defined($al1) && $bo eq $osrc && $bd eq $iaAdId;
        $al2="<Link>$b</Link>" if !defined($al2) && $bo eq $iaAdId && $bd eq $odst;
        last if defined $al1 && defined $al2; }
      die "aproj: adapter links not found for $tmplName\n" unless defined $al1 && defined $al2;
      my $a1=reroll_ids($al1,$ADL[$adi*2]);
      $a1 =~ s/<DestinationID>\Q$iaAdId\E<\/DestinationID>/<DestinationID>$AD[$adi]<\/DestinationID>/ or die;
      my $a2=reroll_ids($al2,$ADL[$adi*2+1]);
      $a2 =~ s/<OriginID>\Q$iaAdId\E<\/OriginID>/<OriginID>$AD[$adi]<\/OriginID>/ or die;
      $a2 =~ s/<DestinationID>\Q$odst\E<\/DestinationID>/<DestinationID>$N{$key}<\/DestinationID>/ or die;
      substr($x,index($x,$al2)+length($al2),0)="\n".$a1."\n".$a2;
      add_to_list(\$x,$osrc,$ADL[$adi*2],'OutLinks');
      add_to_list(\$x,$AD[$adi],$ADL[$adi*2],'InLinks');
      add_to_list(\$x,$AD[$adi],$ADL[$adi*2+1],'OutLinks');
      add_to_list(\$x,$N{$key},$ADL[$adi*2+1],'InLinks'); }
    print "aproj: +8 adapter links\n"; }

  # 24 synapse links
  { my $tmpl = find_link(\$x,$R{tmplE},$R{tmpltD}); die "aproj aff tmpl\n" unless $tmpl;
    my $at=index($x,$tmpl)+length($tmpl); my $ins='';
    for my $i (0..$#AFF) {
      my ($src,$dst)=@{$AFF[$i]};
      my $new=reroll_ids($tmpl,$LK[$i]);
      $new =~ s/<OriginID>\Q$R{tmplE}\E<\/OriginID>/<OriginID>$R{$src}<\/OriginID>/ or die;
      $new =~ s/<DestinationID>\Q$R{tmpltD}\E<\/DestinationID>/<DestinationID>$R{$dst}<\/DestinationID>/ or die;
      $new =~ s/<SynapticTypeID>[^<]*<\/SynapticTypeID>/<SynapticTypeID>$TYPE<\/SynapticTypeID>/ or die;
      $ins .= $new."\n"; }
    substr($x,$at,0)="\n".$ins;
    for my $i (0..$#AFF) {
      my ($src,$dst)=@{$AFF[$i]};
      add_to_list(\$x,$R{$src},$LK[$i],'OutLinks');
      add_to_list(\$x,$R{$dst},$LK[$i],'InLinks'); }
    print "aproj: +",scalar(@AFF)," afferent links\n"; }

  # drawings: II nodes + adapters near their Ia counterparts; links cloned
  { my ($cbeg)=$x =~ /<DiagramXml><!\[CDATA\[/ or die "no CDATA\n";
    my $cs=$cbeg; my $ce=index($x,']]></DiagramXml>',$cs);
    my $cd=substr($x,$cs,$ce-$cs);
    my ($tmplnode)=$cd =~ m{(<Node Left="[^"]+" Top="[^"]+"[^>]*>(?:(?!</Node>).)*?<Tag>\Q$R{iaHipL}\E</Tag>(?:(?!</Node>).)*?</Node>)}s;
    die "II drawing tmpl\n" unless $tmplnode;
    my @names=('L Hip flx II','R Hip flx II','L Knee flx II','R Knee flx II');
    my $ins='';
    for my $i (0..3) { my $n=$tmplnode; my $top=560+$i*50;
      $n =~ s/Left="[-\d.]+"/Left="620"/; $n =~ s/Top="[-\d.]+"/Top="$top"/;
      $n =~ s/<Text>[^<]*<\/Text>/<Text>$names[$i]<\/Text>/;
      $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$N{$II[$i][0]}<\/Tag>/;
      $ins .= "  ".$n."\n"; }
    for my $i (0..3) { my $n=$tmplnode; my $top=560+$i*50;
      $n =~ s/Left="[-\d.]+"/Left="690"/; $n =~ s/Top="[-\d.]+"/Top="$top"/;
      $n =~ s/<Text>[^<]*<\/Text>/<Text>$names[$i] Adapter<\/Text>/;
      $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$AD[$i]<\/Tag>/;
      $ins .= "  ".$n."\n"; }
    # link drawings: reuse the aff template link's drawing per new link
    my ($tib) = find_link(\$x,$R{tmplE},$R{tmpltD}) =~ /<ID>([^<]+)/;
    my ($tl) = $cd =~ m{(<Link Org="\d+" Dst="\d+">(?:(?!</Link>).)*?<Tag>\Q$tib\E</Tag>(?:(?!</Link>).)*?</Link>)}s;
    for my $i (0..$#AFF) { my $n=$tl; $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$LK[$i]<\/Tag>/; $ins .= "  $n\n"; }
    for my $i (0..7) { my $n=$tl; $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$ADL[$i]<\/Tag>/; $ins .= "  $n\n"; }
    $cd =~ s{(\s*</AddFlow>)}{\n$ins$1}s or die;
    my $nn=()=$cd=~/<Node Left=/g; my $ln=()=$cd=~/<Link Org=/g;
    $cd =~ s{<AddFlow Nodes="\d+" Links="\d+"}{<AddFlow Nodes="$nn" Links="$ln"};
    print "aproj page: Nodes=$nn Links=$ln\n";
    substr($x,$cs,$ce-$cs)=$cd; }

  spew($PROJ,$x); print "aproj written\n";
}
print "BUILD AFF DONE\n";
