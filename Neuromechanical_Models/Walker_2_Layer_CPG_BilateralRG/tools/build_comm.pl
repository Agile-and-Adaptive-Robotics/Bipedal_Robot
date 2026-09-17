#!/usr/bin/perl
# build_comm.pl — Walker_2_Layer_CPG_BilateralRG milestone 2: Shinohara-style
# commissural pathways (Shinohara et al. 2025 bioRxiv, Aoi/Rybak/Danner).
#
# Adds 4 relay interneurons (cloned from L RG ext IN physiology):
#   L c1, R c1  — crossed flexor inhibition
#   L V3, R V3  — crossed extension excitation
# and 8 synapses (all G=0.5 except V3 outputs at G=0.15, per Shinohara
# Table A.1 weights 1/1/0.5/0.15 scaled to the W2L G=0.5 convention):
#   L RG flx --RG Excite--> L c1 --RG Inhibit--> R RG flx   (c1 pathway)
#   R RG flx --RG Excite--> R c1 --RG Inhibit--> L RG flx
#   L RG ext --RG Excite--> L V3 --RG Excite(G .15)--> R RG ext IN  (V3 pathway)
#   R RG ext --RG Excite--> R V3 --RG Excite(G .15)--> L RG ext IN
#
# Also adds L/R c1 + V3 columns to the Rhythm Generator chart (+aform).
# Purely additive; run AFTER build_rg.pl. Usage: perl build_comm.pl <dir>
use strict; use warnings;

my $dir = shift @ARGV or die "usage: build_comm.pl <dir>\n";
my $PROJ = "$dir/Walker_2_Layer_CPG_BilateralRG.aproj";
my $ASIM = "$dir/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim";
my $AFRM = "$dir/Rhythm_Generator.aform";

sub slurp { open my $f,'<',$_[0] or die "$_[0]: $!"; local $/; my $s=<$f>; close $f; return $s; }
sub spew  { open my $f,'>',$_[0] or die "$_[0]: $!"; print {$f} $_[1]; close $f; }

sub g { my $n=shift; return sprintf("cafe%04d-0000-4000-8000-%012d",$n,900000000000+$n); }
my $childseq = 700;
sub childguid { return g($childseq++); }
sub reroll_ids {
  my ($blk,$main)=@_;
  my $first = 1;
  $blk =~ s{<ID>([^<]+)</ID>}{ $first ? do { $first=0; "<ID>$main</ID>" } : "<ID>".childguid()."</ID>" }ge;
  return $blk;
}

my %N = (Lc1=>g(5), Rc1=>g(6), LV3=>g(7), RV3=>g(8));
my %P = ( LRGext=>'17fdb1f6', LRGextIN=>'6e188aac', LRGflx=>'7bd02f2b', LRGflxIN=>'136ad874',
          RRGext=>'cafe0001', RRGextIN=>'cafe0002', RRGflx=>'cafe0003', RRGflxIN=>'cafe0004' );
my @LINKS = (   # [newid, tmplSrc, tmplDst, newSrc, newDst, G]
  [g(22),$P{LRGext}, $P{LRGextIN}, 'LRGflx',  'Lc1',     0.5 ],   # flx -> c1 (RG Excite)
  [g(23),$P{LRGextIN},$P{LRGflx},  'Lc1',     'RRGflx',  0.5 ],   # c1 -> contra flx (RG Inhibit)
  [g(24),$P{LRGext}, $P{LRGextIN}, 'RRGflx',  'Rc1',     0.5 ],
  [g(25),$P{LRGextIN},$P{LRGflx},  'Rc1',     'LRGflx',  0.5 ],
  [g(26),$P{LRGext}, $P{LRGextIN}, 'LRGext',  'LV3',     0.5 ],   # ext -> V3 (RG Excite)
  [g(27),$P{LRGext}, $P{LRGextIN}, 'LV3',     'RRGextIN',0.15],   # V3 -> contra IN-E (RG Excite, weak)
  [g(28),$P{LRGext}, $P{LRGextIN}, 'RRGext',  'RV3',     0.5 ],
  [g(29),$P{LRGext}, $P{LRGextIN}, 'RV3',     'LRGextIN',0.15],
);
my @NEUR = ( ['Lc1','L c1'], ['Rc1','R c1'], ['LV3','L V3'], ['RV3','R V3'] );
sub resolve {
  my ($xml,$pfx) = @_;
  my $i = index($xml,"<ID>$pfx");
  die "prefix $pfx not found\n" if $i<0;
  return substr($xml,$i+4,index($xml,'</ID>',$i+4)-$i-4);
}

# ============================ ASIM ============================
{
  my $x = slurp($ASIM);
  die "asim already has commissural layer\n" if $x=~/<Name>L c1<\/Name>/;
  my %R = ( (map { $_ => resolve($x,$P{$_}) } keys %P), %N );

  for my $spec (@NEUR) {
    my ($key,$name)=@$spec;
    my $tid = resolve($x,$P{LRGextIN});          # physiology template
    my $idx = index($x,"<ID>$tid<\/ID>");
    die "asim: template IN neuron not found\n" if $idx<0;
    my $s = rindex($x,'<Neuron>',$idx); my $e = index($x,'</Neuron>',$idx)+9;
    my $blk = substr($x,$s,$e-$s);
    my $new = reroll_ids($blk, $N{$key});
    $new =~ s/<Name>[^<]*<\/Name>/<Name>$name<\/Name>/ or die;
    substr($x,$e,0) = "\n".$new;
    print "asim: +neuron $name\n";
  }

  for my $L (@LINKS) {
    my ($newid,$ts,$td,$ns,$nd,$gg)=@$L;
    my ($to,$tdf) = (resolve($x,$ts), resolve($x,$td));
    my ($found); pos($x)=0;
    while ($x =~ /<Connexion>(.*?)<\/Connexion>/gs) {
      my $b=$1; my ($o)=$b=~/<SourceID>([^<]+)/; my ($d)=$b=~/<TargetID>([^<]+)/;
      if ($o eq $to && $d eq $tdf) { $found=$b; last; }
    }
    die "asim: template connexion $ts->$td not found\n" unless $found;
    my $blk = "<Connexion>$found<\/Connexion>";
    my $new = reroll_ids($blk, $newid);
    $new =~ s/<SourceID>\Q$to\E<\/SourceID>/<SourceID>$R{$ns}<\/SourceID>/ or die;
    $new =~ s/<TargetID>\Q$tdf\E<\/TargetID>/<TargetID>$R{$nd}<\/TargetID>/ or die;
    $new =~ s/<G>[^<]*<\/G>/<G>$gg<\/G>/ or die "asim: G patch failed\n";
    substr($x,index($x,$blk)+length($blk),0) = "\n".$new;
    print "asim: +connexion $ns -> $nd G=$gg\n";
  }

  # chart columns
  {
    my $at = index($x,"<OutputFilename>Rhythm Generator.txt<\/OutputFilename>");
    die "asim: RG chart missing\n" if $at<0;
    my $dcend = index($x,"<\/DataColumns>",$at);
    my $t = index($x,"<ColumnName>R RG flx<\/ColumnName>");
    die "asim: R RG flx column missing (run build_rg first)\n" if $t<0;
    my $ts2 = rindex($x,'<DataColumn>',$t); my $te = index($x,'</DataColumn>',$t)+13;
    my $tmplblk = substr($x,$ts2,$te-$ts2);
    my @col = ( [g(35),'L c1',$N{Lc1}], [g(36),'R c1',$N{Rc1}], [g(37),'L V3',$N{LV3}], [g(38),'R V3',$N{RV3}] );
    my $ins='';
    for my $c (@col) {
      my ($cid,$cn,$tgt)=@$c;
      my $new=$tmplblk;
      my ($oldid)=$tmplblk=~/<ID>([^<]+)/;
      $new =~ s/<ID>\Q$oldid\E<\/ID>/<ID>$cid<\/ID>/ or die;
      $new =~ s/<ColumnName>[^<]+<\/ColumnName>/<ColumnName>$cn<\/ColumnName>/ or die;
      $new =~ s/<TargetID>[^<]+<\/TargetID>/<TargetID>$tgt<\/TargetID>/ or die;
      $ins .= $new."\n";
    }
    substr($x,$dcend,0) = $ins;
    print "asim: +4 commissural chart columns\n";
  }
  spew($ASIM,$x);
  print "asim written.\n";
}

# ============================ APROJ ============================
{
  my $x = slurp($PROJ);
  die "aproj already has commissural layer\n" if $x=~/<Text>L c1<\/Text>/;
  my %R = ( (map { $_ => resolve($x,$P{$_}) } keys %P), %N );

  sub node_block { my ($x,$id)=@_;
    my $idx = index($x,"<ID>$id<\/ID>"); return () if $idx<0;
    my $s = rindex($x,'<Node>',$idx); my $e = index($x,'</Node>',$idx)+7;
    return () if $s<0 || $e<7; return ($s,$e-$s); }
  sub add_to_list { my ($xr,$nid,$lid,$list)=@_;
    my ($s,$l)=node_block($$xr,$nid); die "node $nid not found for list add\n" unless $s;
    my $seg = substr($$xr,$s,$l); my $nseg = $seg;
    my $rep = "<$list>\n<ID>$lid<\/ID>\n</$list>";
    if    ($nseg =~ s{<$list/>\s*}{$rep}s) {}
    elsif ($nseg =~ s{<$list>\s*</$list>}{$rep}s) {}
    elsif ($nseg =~ s{</$list>}{<ID>$lid<\/ID>\n</$list>}s) {}
    else { die "no $list on node $nid\n"; }
    substr($$xr,$s,$l)=$nseg; }
  sub find_synlink { my ($x,$o,$d)=@_; pos($$x)=0;
    while ($$x =~ /<Link>(.*?)<\/Link>/gs) {
      my $b=$1; next unless $b =~ /Behavior\.Synapse</;
      my ($bo)=$b=~/<OriginID>([^<]+)/; my ($bd)=$b=~/<DestinationID>([^<]+)/;
      return "<Link>$b<\/Link>" if $bo eq $o && $bd eq $d; }
    return undef; }

  # new nodes
  for my $spec (@NEUR) {
    my ($key,$name)=@$spec;
    my $tid = resolve($x,$P{LRGextIN});
    my ($s,$l)=node_block($x,$tid); die "aproj: template IN node missing\n" unless $s;
    my $blk = substr($x,$s,$l);
    my $new = reroll_ids($blk, $N{$key});
    $new =~ s/<Text>[^<]*<\/Text>/<Text>$name<\/Text>/ or die;
    $new =~ s{<InLinks>.*?</InLinks>}{<InLinks>\n</InLinks>}gs;
    $new =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks>\n</OutLinks>}gs;
    substr($x,$s+$l,0) = "\n".$new;
    print "aproj: +node $name\n";
  }

  # new links
  for my $L (@LINKS) {
    my ($newid,$ts,$td,$ns,$nd,$gg)=@$L;
    my ($to,$tdf)=(resolve($x,$ts),resolve($x,$td));
    my $blk = find_synlink(\$x,$to,$tdf);
    die "aproj: template synapse $ts->$td not found\n" unless $blk;
    my $new = reroll_ids($blk, $newid);
    $new =~ s/<OriginID>\Q$to\E<\/OriginID>/<OriginID>$R{$ns}<\/OriginID>/ or die;
    $new =~ s/<DestinationID>\Q$tdf\E<\/DestinationID>/<DestinationID>$R{$nd}<\/DestinationID>/ or die;
    # G patch: aproj carries conductance in SynapticConductance Value/Actual (micro S)
    my $actual = sprintf("%.3g", $gg*1e-6); $actual =~ s/e-0/e-/;
    $new =~ s{<SynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>}{<SynapticConductance Value="$gg" Scale="micro" Actual="$actual"/>}
      or die "aproj: conductance patch failed\n";
    substr($x,index($x,$blk)+length($blk),0) = "\n".$new;
    add_to_list(\$x,$R{$ns},$newid,'OutLinks');
    add_to_list(\$x,$R{$nd},$newid,'InLinks');
    print "aproj: +link $ns -> $nd G=$gg\n";
  }

  # drawings
  {
    my ($cbeg) = $x =~ /<DiagramXml><!\[CDATA\[/ or die "no page CDATA\n";
    my $cs=$cbeg; my $ce = index($x,']]></DiagramXml>',$cs);
    my $cd = substr($x,$cs,$ce-$cs);
    my $tmplid = $R{LRGext};
    my ($tmplnode) = $cd =~ m{(<Node Left="[^"]+" Top="[^"]+"[^>]*>(?:(?!</Node>).)*?<Tag>\Q$tmplid\E</Tag>(?:(?!</Node>).)*?</Node>)}s
      or die "L RG ext drawing entry not found\n";
    my @pos = ( [100,900],[280,900],[180,960],[360,960] );
    my @names = ('L c1','R c1','L V3','R V3');
    my @keys = ('Lc1','Rc1','LV3','RV3');
    my $ins='';
    for my $i (0..3) {
      my $n=$tmplnode;
      $n =~ s/Left="[-\d.]+"/Left="$pos[$i][0]"/;
      $n =~ s/Top="[-\d.]+"/Top="$pos[$i][1]"/;
      $n =~ s/<Text>[^<]*<\/Text>/<Text>$names[$i]<\/Text>/;
      $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$N{$keys[$i]}<\/Tag>/;
      $ins .= "  ".$n."\n";
    }
    $cd =~ s{(\s*</AddFlow>)}{\n$ins$1}s or die "AddFlow close missing\n";
    for my $L (@LINKS) {
      my ($newid,$ts,$td,$ns,$nd,$gg)=@$L;
      my ($to,$tdf)=(resolve($x,$ts),resolve($x,$td));
      my $fblk = find_synlink(\$x,$to,$tdf);
      # NB: find_synlink finds the ORIGINAL template (first match) — but clones
      # of the same template pair exist now. Match by (o,d) returns the first
      # (original L-side) block, whose drawing entry is the right template.
      die "template vanished\n" unless $fblk;
      my ($tib)=$fblk=~/<ID>([^<]+)/;
      my ($tl) = $cd =~ m{(<Link Org="\d+" Dst="\d+">(?:(?!</Link>).)*?<Tag>\Q$tib\E</Tag>(?:(?!</Link>).)*?</Link>)}s;
      unless ($tl) { print "WARN: no drawing for template $tib; skip drawing\n"; next; }
      my $dn=$tl; $dn =~ s/<Tag>[^<]*<\/Tag>/<Tag>$newid<\/Tag>/;
      $cd =~ s{(\s*</AddFlow>)}{\n  $dn$1}s;
    }
    my $nn=()=$cd=~/<Node Left=/g; my $ln=()=$cd=~/<Link Org=/g;
    $cd =~ s{<AddFlow Nodes="\d+" Links="\d+"}{<AddFlow Nodes="$nn" Links="$ln"} or die;
    print "aproj: page now Nodes=$nn Links=$ln\n";
    substr($x,$cs,$ce-$cs)=$cd;
  }
  spew($PROJ,$x);
  print "aproj written.\n";
}

# ============================ AFORM ============================
{
  my $x = slurp($AFRM);
  die "aform already has commissural columns\n" if $x=~/<Name>L c1<\/Name>/;
  my ($t) = $x =~ /(<DataColumn>(?:(?!<\/DataColumn>).)*?<Name>R RG flx<\/Name>(?:(?!<\/DataColumn>).)*?<\/DataColumn>)/s
    or die "aform: R RG flx column not found\n";
  my $tpos = index($x,"<Name>R RG flx<\/Name>");
  my $dcend = index($x,"<\/DataColumns>",$tpos);
  die "aform: DataColumns close missing\n" if $dcend<0;
  my @col = ( [g(35),'L c1',$N{Lc1},-65536], [g(36),'R c1',$N{Rc1},-16711936],
              [g(37),'L V3',$N{LV3},-65281], [g(38),'R V3',$N{RV3},-256] );
  my $ins='';
  for my $c (@col) {
    my ($cid,$cn,$tgt,$col)=@$c;
    my $new=$t;
    my ($oldid)=$t=~/<ID>([^<]+)/;
    $new =~ s/<ID>\Q$oldid\E<\/ID>/<ID>$cid<\/ID>/ or die;
    $new =~ s/<Name>[^<]*<\/Name>/<Name>$cn<\/Name>/ or die;
    $new =~ s/<DataItemID>[^<]+<\/DataItemID>/<DataItemID>$tgt<\/DataItemID>/ or die;
    $new =~ s/<LineColor>-?\d+<\/LineColor>/<LineColor>$col<\/LineColor>/;
    $ins .= $new."\n";
  }
  substr($x,$dcend,0) = "\n".$ins;
  spew($AFRM,$x);
  print "aform written (+4 columns).\n";
}
print "BUILD COMM DONE\n";
