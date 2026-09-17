#!/usr/bin/perl
# build_rg.pl — Walker_2_Layer_CPG_BilateralRG milestone 1: bilateral RG.
#
# Transforms (in lockstep, shared deterministic GUIDs, prefix cafe____):
#   1. Remove the 4 crossed L RG -> R PF connexions (single-RG antiphase wiring).
#   2. Add R RG half-center: R RG ext / R RG ext IN / R RG flx / R RG flx IN
#      (cloned from L counterparts), wired ipsilaterally:
#        ext -> ext IN (RG Excite), ext IN -> flx (RG Inhibit),
#        flx -> flx IN (RG Excite), flx IN -> ext (RG Inhibit),
#        ext -> flx (RG to RG Excite), flx -> ext (RG to RG Excite),
#        ext -> R Hip PF ext, ext -> R Knee PF ext (RG to PF Excite),
#        flx -> R Hip PF flx, flx -> R Knee PF flx (RG to PF Excite).   All G=0.5.
#   3. Stimulus_2: 10 nA tonic current 0..0.01 s into R RG flx
#      (Stimulus_1 already kicks L RG ext -> legs start neurally antiphase).
#   4. Rhythm Generator chart (+aform): 4 R RG membrane-voltage columns.
#   5. Start pose: femur_L/R Rotation Z = hipL/hipR deg, tibia_L/R Rotation Y
#      = kneeL/kneeR deg (degrees in aproj, radians in asim).
#
# Files touched: <dir>/Walker_2_Layer_CPG_BilateralRG.aproj
#                <dir>/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim
#                <dir>/Rhythm_Generator.aform
# Usage: perl build_rg.pl <dir> [hipL hipR kneeL kneeR]   (deg, default 12 -12 -28 -28)
# Refuses to run twice (checks for existing "R RG ext" neuron).
use strict; use warnings;

my $dir = shift @ARGV or die "usage: build_rg.pl <dir> [hipL hipR kneeL kneeR]\n";
my ($hipL,$hipR,$kneeL,$kneeR) = (@ARGV, 12, -12, -28, -28)[0..3];
my $PROJ = "$dir/Walker_2_Layer_CPG_BilateralRG.aproj";
my $ASIM = "$dir/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim";
my $AFRM = "$dir/Rhythm_Generator.aform";

sub slurp { open my $f,'<',$_[0] or die "$_[0]: $!"; local $/; my $s=<$f>; close $f; return $s; }
sub spew  { open my $f,'>',$_[0] or die "$_[0]: $!"; print {$f} $_[1]; close $f; }

# ---- deterministic GUIDs ----
sub g { my $n=shift; return sprintf("cafe%04d-0000-4000-8000-%012d",$n,900000000000+$n); }
my $childseq = 400;
sub childguid { return g($childseq++); }
# Replace the FIRST <ID> in a cloned block with $main, re-roll every other
# child-object <ID> (Ca channels, gains, ...) so clones stay unique.
sub reroll_ids {
  my ($blk,$main)=@_;
  my $first = 1;
  $blk =~ s{<ID>([^<]+)</ID>}{ $first ? do { $first=0; "<ID>$main</ID>" } : "<ID>".childguid()."</ID>" }ge;
  return $blk;
}
my %N = (RRGext=>g(1), RRGextIN=>g(2), RRGflx=>g(3), RRGflxIN=>g(4));
my @LK = map { g(11+$_) } 0..9;          # 10 new synapse links
my $STIM2 = g(21);
my @CL = map { g(31+$_) } 0..3;          # 4 chart columns

# ---- known existing IDs (8-char prefixes, resolved per file) ----
my %P = (
  LRGext=>'17fdb1f6', LRGextIN=>'6e188aac', LRGflx=>'7bd02f2b', LRGflxIN=>'136ad874',
  LPFhipE=>'57b344bf', LPFhipF=>'f5ba89bb', LPFkneeE=>'f336ecdc', LPFkneeF=>'9bc35128',
  RPFhipE=>'aa4d0055', RPFhipF=>'16bc9942', RPFkneeE=>'eb127c29', RPFkneeF=>'fbda99c8',
);
sub resolve {
  my ($xml,$pfx) = @_;
  my $i = index($xml,"<ID>$pfx");
  die "prefix $pfx not found\n" if $i<0;
  return substr($xml,$i+4,index($xml,'</ID>',$i+4)-$i-4);
}

# ---- link spec: [newguid, tmplSrc, tmplDst, newSrcKey, newDstKey] ----
my @LINKS = (
  [$LK[0],$P{LRGext},  $P{LRGextIN}, 'RRGext',  'RRGextIN'],
  [$LK[1],$P{LRGextIN},$P{LRGflx},   'RRGextIN','RRGflx'],
  [$LK[2],$P{LRGflx},  $P{LRGflxIN}, 'RRGflx',  'RRGflxIN'],
  [$LK[3],$P{LRGflxIN},$P{LRGext},   'RRGflxIN','RRGext'],
  [$LK[4],$P{LRGext},  $P{LRGflx},   'RRGext',  'RRGflx'],   # RG to RG
  [$LK[5],$P{LRGflx},  $P{LRGext},   'RRGflx',  'RRGext'],   # RG to RG
  [$LK[6],$P{LRGext},  $P{LPFhipE},  'RRGext',  'RPFhipE'],
  [$LK[7],$P{LRGext},  $P{LPFkneeE}, 'RRGext',  'RPFkneeE'],
  [$LK[8],$P{LRGflx},  $P{LPFhipF},  'RRGflx',  'RPFhipF'],
  [$LK[9],$P{LRGflx},  $P{LPFkneeF}, 'RRGflx',  'RPFkneeF'],
);
# ---- crossed links to delete: [src,dst] ----
my @CROSS = ( [$P{LRGext},$P{RPFhipF}], [$P{LRGext},$P{RPFkneeF}],
              [$P{LRGflx},$P{RPFhipE}], [$P{LRGflx},$P{RPFkneeE}] );

my @NEUR = ( ['RRGext',  $P{LRGext},  'R RG ext'],
             ['RRGextIN',$P{LRGextIN},'R RG ext IN'],
             ['RRGflx',  $P{LRGflx},  'R RG flx'],
             ['RRGflxIN',$P{LRGflxIN},'R RG flx IN'] );

my $pi = 4*atan2(1,1);
sub d2r { sprintf("%.7f", $_[0]*$pi/180) }

# ============================ ASIM ============================
{
  my $x = slurp($ASIM);
  die "asim already has R RG ext (refusing double-run)\n" if $x=~/<Name>R RG ext<\/Name>/;
  my %R = ( (map { $_ => resolve($x,$P{$_}) } keys %P), %N );

  # --- neurons ---
  for my $spec (@NEUR) {
    my ($key,$tmpl,$name)=@$spec;
    my $tid = resolve($x,$tmpl);
    my $idx = index($x,"<ID>$tid<\/ID>");
    die "asim: template neuron $tmpl not found\n" if $idx<0;
    my $s = rindex($x,'<Neuron>',$idx); my $e = index($x,'</Neuron>',$idx)+9;
    die "asim: neuron block bounds failed for $tmpl\n" if $s<0||$e<9;
    my $blk = substr($x,$s,$e-$s);
    my $new = reroll_ids($blk, $N{$key});
    $new =~ s/<Name>[^<]*<\/Name>/<Name>$name<\/Name>/ or die;
    substr($x,$e,0) = "\n".$new;
    print "asim: +neuron $name\n";
  }

  # --- connexions: clone 10 ---
  for my $L (@LINKS) {
    my ($newid,$ts,$td,$ns,$nd)=@$L;
    my ($to,$tdf) = (resolve($x,$ts), resolve($x,$td));
    my ($found);
    pos($x)=0;
    while ($x =~ /<Connexion>(.*?)<\/Connexion>/gs) {
      my $b=$1; my ($o)=$b=~/<SourceID>([^<]+)/; my ($d)=$b=~/<TargetID>([^<]+)/;
      if ($o eq $to && $d eq $tdf) { $found=$b; last; }
    }
    die "asim: template connexion $ts->$td not found\n" unless $found;
    my $blk = "<Connexion>$found<\/Connexion>";
    my $new = reroll_ids($blk, $newid);
    $new =~ s/<SourceID>\Q$to\E<\/SourceID>/<SourceID>$R{$ns}<\/SourceID>/ or die;
    $new =~ s/<TargetID>\Q$tdf\E<\/TargetID>/<TargetID>$R{$nd}<\/TargetID>/ or die;
    my $at = index($x,$blk)+length($blk);
    substr($x,$at,0) = "\n".$new;
    print "asim: +connexion $ns -> $nd ($newid)\n";
  }

  # --- connexions: delete 4 crossed ---
  for my $C (@CROSS) {
    my ($s,$d)=(resolve($x,$C->[0]),resolve($x,$C->[1]));
    my $done=0;
    pos($x)=0;
    while ($x =~ /<Connexion>(.*?)<\/Connexion>/gs) {
      my $b=$1; my ($o)=$b=~/<SourceID>([^<]+)/; my ($d2)=$b=~/<TargetID>([^<]+)/;
      next unless $o eq $s && $d2 eq $d;
      my $blk="<Connexion>$b<\/Connexion>";
      my $at=index($x,$blk);
      my $del=$blk."\n";
      substr($x,$at,length($del))="";
      print "asim: -crossed connexion $C->[0] -> $C->[1]\n"; $done=1; last;
    }
    die "asim: crossed link $C->[0]->$C->[1] not found\n" unless $done;
  }

  # --- Stimulus_2 ---
  {
    my $at = index($x,"<Name>Stimulus_1<\/Name>");
    die "asim: Stimulus_1 not found\n" if $at<0;
    my $s = rindex($x,'<Stimulus>',$at); my $e = index($x,'</Stimulus>',$at)+11;
    my $blk = substr($x,$s,$e-$s);
    my $new = reroll_ids($blk, $STIM2);
    $new =~ s/<Name>Stimulus_1<\/Name>/<Name>Stimulus_2<\/Name>/ or die;
    $new =~ s/<TargetNodeID>[^<]+<\/TargetNodeID>/<TargetNodeID>$N{RRGflx}<\/TargetNodeID>/ or die;
    substr($x,$e,0) = "\n".$new;
    print "asim: +Stimulus_2 -> R RG flx\n";
  }

  # --- chart columns (Rhythm Generator DataChart) ---
  {
    my $at = index($x,"<OutputFilename>Rhythm Generator.txt<\/OutputFilename>");
    die "asim: RG chart not found\n" if $at<0;
    my $dcend = index($x,"<\/DataColumns>",$at);
    die "asim: RG DataColumns end not found\n" if $dcend<0;
    # template: the L RG ext column
    my $t = index($x,"<ColumnName>L RG ext<\/ColumnName>");
    die "asim: L RG ext column not found\n" if $t<0;
    my $ts = rindex($x,'<DataColumn>',$t); my $te = index($x,'</DataColumn>',$t)+13;
    my $tmplblk = substr($x,$ts,$te-$ts);
    my @colspec = ( [$CL[0],'R RG ext IN',$N{RRGextIN}], [$CL[1],'R RG flx IN',$N{RRGflxIN}],
                    [$CL[2],'R RG ext',$N{RRGext}],      [$CL[3],'R RG flx',$N{RRGflx}] );
    my $ins='';
    for my $c (@colspec) {
      my ($cid,$cn,$tgt)=@$c;
      my $new=$tmplblk;
      my ($oldid)=$tmplblk=~/<ID>([^<]+)/;
      $new =~ s/<ID>\Q$oldid\E<\/ID>/<ID>$cid<\/ID>/ or die;
      $new =~ s/<ColumnName>[^<]+<\/ColumnName>/<ColumnName>$cn<\/ColumnName>/ or die;
      $new =~ s/<TargetID>[^<]+<\/TargetID>/<TargetID>$tgt<\/TargetID>/ or die;
      $ins .= $new."\n";
    }
    substr($x,$dcend,0) = $ins;
    print "asim: +4 RG chart columns\n";
  }

  # --- start pose (radians in asim) ---
  #  femur_L/R: <Rotation x y z/> z=hip ; tibia_L/R: y=knee
  my @pose = ( ['femur_L','z',d2r($hipL)], ['femur_R','z',d2r($hipR)],
               ['tibia_L','y',d2r($kneeL)], ['tibia_R','y',d2r($kneeR)] );
  for my $p (@pose) {
    my ($body,$ax,$val)=@$p;
    my $at = index($x,"<Name>$body<\/Name>");
    die "asim: body $body not found\n" if $at<0;
    my $bs = rindex($x,'<RigidBody>',$at); my $be = index($x,'</RigidBody>',$at);
    my $seg = substr($x,$bs,$be-$bs);
    my ($rline) = $seg =~ /(<Rotation [^>]+\/>)/ or die "asim: $body rotation line not found\n";
    my $nl = $rline;
    if    ($ax eq 'z') { $nl =~ s/z="[-\d.e]+"/z="$val"/ or die; }
    else               { $nl =~ s/y="[-\d.e]+"/y="$val"/ or die; }
    substr($x,$bs,$be-$bs) =~ s/\Q$rline\E/$nl/ or die "asim: pose replace failed $body\n";
    print "asim: pose $body $ax=$val rad\n";
  }

  spew($ASIM,$x);
  print "asim written.\n";
}

# ============================ APROJ ============================
{
  my $x = slurp($PROJ);
  die "aproj already has R RG ext (refusing double-run)\n" if $x=~/<Text>R RG ext<\/Text>/;
  my %R = ( (map { $_ => resolve($x,$P{$_}) } keys %P), %N );

  # helpers -------------------------------------------------------
  sub node_block {   # returns (start,end) of <Node> block owning this ID
    my ($x,$id)=@_;
    my $idx = index($x,"<ID>$id<\/ID>");
    return () if $idx<0;
    my $s = rindex($x,'<Node>',$idx); my $e = index($x,'</Node>',$idx)+7;
    return () if $s<0 || $e<7;
    return ($s,$e-$s);
  }
  sub add_to_list {  # append link id to node's InLinks/OutLinks
    my ($xr,$nid,$lid,$list)=@_;
    my ($s,$l)=node_block($$xr,$nid); die "node $nid not found for list add\n" unless $s;
    my $seg = substr($$xr,$s,$l);
    my $nseg = $seg;
    my $rep = "<$list>\n<ID>$lid<\/ID>\n</$list>";
    if    ($nseg =~ s{<$list/>\s*}{$rep}s) {}
    elsif ($nseg =~ s{<$list>\s*</$list>}{$rep}s) {}
    elsif ($nseg =~ s{</$list>}{<ID>$lid<\/ID>\n</$list>}s) {}
    else { die "no $list on node $nid\n"; }
    substr($$xr,$s,$l)=$nseg;
  }
  sub strip_from_list {
    my ($xr,$nid,$lid)=@_;
    my ($s,$l)=node_block($$xr,$nid); die "node $nid not found for list strip\n" unless $s;
    my $seg = substr($$xr,$s,$l);
    my $nseg = $seg;
    $nseg =~ s{<ID>\Q$lid\E<\/ID>\r?\n?}{}g;
    substr($$xr,$s,$l)=$nseg;
  }
  my $nlcount=0;
  sub find_synlink {  # functional <Link> block by origin/dest (Behavior.Synapse)
    my ($x,$o,$d)=@_;
    pos($$x)=0;
    while ($$x =~ /<Link>(.*?)<\/Link>/gs) {
      my $b=$1;
      next unless $b =~ /Behavior\.Synapse</;
      my ($bo)=$b=~/<OriginID>([^<]+)/; my ($bd)=$b=~/<DestinationID>([^<]+)/;
      return "<Link>$b<\/Link>" if $bo eq $o && $bd eq $d;
    }
    return undef;
  }

  # --- 1. delete crossed links first (so lists clean before adds) ---
  my @crossIds;
  for my $C (@CROSS) {
    my ($s,$d)=(resolve($x,$C->[0]),resolve($x,$C->[1]));
    my $blk = find_synlink(\$x,$s,$d);
    die "aproj: crossed link $C->[0]->$C->[1] not found\n" unless $blk;
    my ($lid)=$blk=~/<ID>([^<]+)/;
    push @crossIds,$lid;
    substr($x,index($x,$blk),length($blk))='';
    strip_from_list(\$x,$s,$lid);
    strip_from_list(\$x,$d,$lid);
    print "aproj: -crossed link $C->[0] -> $C->[1] ($lid)\n";
  }

  # --- 2. new nodes ---
  for my $spec (@NEUR) {
    my ($key,$tmpl,$name)=@$spec;
    my $tid = resolve($x,$tmpl);
    my ($s,$l)=node_block($x,$tid); die "aproj: template node $tmpl missing\n" unless $s;
    my $blk = substr($x,$s,$l);
    my $new = reroll_ids($blk, $N{$key});
    $new =~ s/<Text>[^<]*<\/Text>/<Text>$name<\/Text>/ or die;
    $new =~ s{<InLinks>.*?</InLinks>}{<InLinks>\n</InLinks>}gs;
    $new =~ s{<OutLinks>.*?</OutLinks>}{<OutLinks>\n</OutLinks>}gs;
    substr($x,$s+$l,0) = "\n".$new;
    print "aproj: +node $name\n";
  }

  # --- 3. new synapse links ---
  for my $L (@LINKS) {
    my ($newid,$ts,$td,$ns,$nd)=@$L;
    my ($to,$tdf)=(resolve($x,$ts),resolve($x,$td));
    my $blk = find_synlink(\$x,$to,$tdf);
    die "aproj: template synapse $ts->$td not found\n" unless $blk;
    my $new = reroll_ids($blk, $newid);
    $new =~ s/<OriginID>\Q$to\E<\/OriginID>/<OriginID>$R{$ns}<\/OriginID>/ or die;
    $new =~ s/<DestinationID>\Q$tdf\E<\/DestinationID>/<DestinationID>$R{$nd}<\/DestinationID>/ or die;
    substr($x,index($x,$blk)+length($blk),0) = "\n".$new;
    add_to_list(\$x,$R{$ns},$newid,'OutLinks');
    add_to_list(\$x,$R{$nd},$newid,'InLinks');
    print "aproj: +link $ns -> $nd\n";
  }

  # --- 4. Stimulus_2 ---
  {
    my $at = index($x,"<Name>Stimulus_1<\/Name>");
    die "aproj: Stimulus_1 not found\n" if $at<0;
    my $s = rindex($x,'<Stimulus>',$at); my $e = index($x,'</Stimulus>',$at)+11;
    my $blk = substr($x,$s,$e-$s);
    my $new = reroll_ids($blk, $STIM2);
    $new =~ s/<Name>Stimulus_1<\/Name>/<Name>Stimulus_2<\/Name>/ or die;
    $new =~ s/<NodeID>[^<]+<\/NodeID>/<NodeID>$N{RRGflx}<\/NodeID>/ or die;
    substr($x,$e,0)="\n".$new;
    print "aproj: +Stimulus_2\n";
  }

  # --- 5. start pose (degrees, attribute triplets) ---
  my @pose = ( ['femur_L','Z',$hipL], ['femur_R','Z',$hipR],
               ['tibia_L','Y',$kneeL], ['tibia_R','Y',$kneeR] );
  for my $p (@pose) {
    my ($body,$ax,$val)=@$p;
    my $at = index($x,"<Name>$body<\/Name>");
    die "aproj: body $body not found\n" if $at<0;
    my $bs = rindex($x,'<RigidBody>',$at); my $be = index($x,'</RigidBody>',$at);
    my $seg = substr($x,$bs,$be-$bs);
    # first <Rotation> block directly under the body (before any <Joint>)
    my $j = index($seg,'<Joint>');
    my $head = $j>0 ? substr($seg,0,$j) : $seg;
    $head =~ /<Rotation>(.*?)<\/Rotation>/s or die "aproj: $body Rotation not found\n";
    my $rblk = $1;
    my $nrblk = $rblk;
    $nrblk =~ s{<$ax Value="[^"]*" Scale="([^"]*)" Actual="[^"]*"/>}{<$ax Value="$val" Scale="$1" Actual="$val"/>} or die "aproj: $body $ax pose replace failed\n";
    my $nseg = $seg; $nseg =~ s/\Q$rblk\E/$nrblk/ or die;
    substr($x,$bs,$be-$bs)=$nseg;
    print "aproj: pose $body $ax=$val deg\n";
  }

  # --- 6. drawings on the single page ---
  {
    my ($cbeg) = $x =~ /<DiagramXml><!\[CDATA\[/ or die "no page CDATA\n";
    my $cs = $cbeg; my $ce = index($x,']]></DiagramXml>',$cs);
    my $cd = substr($x,$cs,$ce-$cs);

    # NOTE: crossed drawing links are deliberately KEPT as orphans. Removing
    # drawing <Link> entries corrupts other links' Org/Dst handle resolution
    # (GUI: "Unable to cast MyLink to Lassalle.Flow.Node", verified 2026-09-16).
    # The stale arrows get rebuilt in the subnetwork page reorganization.
    # 6b. append 4 node drawings cloned from L RG ext drawing
    my ($tmplnode) = $cd =~ /(<Node Left="[^"]+" Top="[^"]+"[^>]*>(?:(?!<\/Node>).)*?<Tag>\Q$R{LRGext}\E<\/Tag>(?:(?!<\/Node>).)*?<\/Node>)/s
      or die "L RG ext drawing entry not found\n";
    my @pos = ( [300,760], [220,760], [300,820], [220,820] );  # ext, extIN, flx, flxIN
    my @names = ('R RG ext','R RG ext IN','R RG flx','R RG flx IN');
    my @keys  = ('RRGext','RRGextIN','RRGflx','RRGflxIN');
    my $ins='';
    for my $i (0..3) {
      my $n = $tmplnode;
      $n =~ s/Left="[-\d.]+"/Left="$pos[$i][0]"/;
      $n =~ s/Top="[-\d.]+"/Top="$pos[$i][1]"/;
      $n =~ s/<Text>[^<]*<\/Text>/<Text>$names[$i]<\/Text>/;
      $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$N{$keys[$i]}<\/Tag>/;
      $ins .= "  ".$n."\n";
    }
    $cd =~ s{(\s*</AddFlow>)}{\n$ins$1}s or die "AddFlow close not found\n";
    # 6c. append 10 link drawings: clone each template's drawing (same template link as functional clone), swap Tag
    for my $L (@LINKS) {
      my ($newid,$ts,$td,$ns,$nd)=@$L;
      # find the ORIGINAL template functional link id (L side) to copy its drawing
      my ($to,$tdf)=(resolve($x,$ts),resolve($x,$td));
      # NB: crossed deletion already removed some; find template among remaining by re-resolving? templates are L-side, untouched.
      my $fblk = find_synlink(\$x,$to,$tdf);
      die "template functional link vanished ($ts->$td)\n" unless $fblk;
      my ($tid2)=$fblk=~/<ID>([^<]+)/;
      my ($tmplink) = $cd =~ m{(<Link Org="\d+" Dst="\d+">(?:(?!</Link>).)*?<Tag>\Q$tid2\E</Tag>(?:(?!</Link>).)*?</Link>)}s;
      unless ($tmplink) { print "WARN: no drawing for template $tid2 ($ts->$td); skipping drawing\n"; next; }
      my $n = $tmplink;
      $n =~ s/<Tag>[^<]*<\/Tag>/<Tag>$newid<\/Tag>/;
      $cd =~ s{(\s*</AddFlow>)}{\n  $n$1}s;
    }
    # 6d. counters
    my $nn = () = $cd =~ /<Node\s+Left=/g;
    my $ln = () = $cd =~ /<Link\s+Org=/g;
    $cd =~ s{<AddFlow Nodes="\d+" Links="\d+"}{<AddFlow Nodes="$nn" Links="$ln"} or die "AddFlow counter line not found\n";
    print "aproj: page now Nodes=$nn Links=$ln\n";

    substr($x,$cs,$ce-$cs) = $cd;
  }

  spew($PROJ,$x);
  print "aproj written.\n";
}

# ============================ AFORM ============================
{
  my $x = slurp($AFRM);
  die "aform already has R RG columns\n" if $x=~/<Name>R RG ext<\/Name>/;
  my ($t) = $x =~ /(<DataColumn>(?:(?!<\/DataColumn>).)*?<Name>L RG ext<\/Name>(?:(?!<\/DataColumn>).)*?<\/DataColumn>)/s
    or die "aform: L RG ext column not found\n";
  # insert just before the </DataColumns> that encloses the L RG columns
  my $tpos = index($x,"<Name>L RG ext<\/Name>");
  my $dcend = index($x,"<\/DataColumns>",$tpos);
  die "aform: DataColumns close after L RG ext not found\n" if $dcend<0;
  my @colspec = ( [$CL[0],'R RG ext IN',$N{RRGextIN},-65536], [$CL[1],'R RG flx IN',$N{RRGflxIN},-16711936],
                  [$CL[2],'R RG ext',$N{RRGext},-65281],     [$CL[3],'R RG flx',$N{RRGflx},-256] );
  my $ins='';
  for my $c (@colspec) {
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

print "BUILD RG DONE\n";
