#!/usr/bin/perl
# patch_comm_types.pl — give the commissural outputs dedicated SynapseTypes.
#
# Finding (2026-09-16, air tests): AnimatSimulator's effective conductance for
# these NonSpiking synapses follows the SynapseType SynAmp; the per-connexion
# <G> in the asim is inert (G=1e-4 vs 0.15 gave bit-identical runs). With
# "RG Excite" (SynAmp 2.749) the V3->contra-IN-E path tonically suppresses
# the contralateral flexor and latches both RGs. With SynAmp 0.1 the network
# oscillates in stable antiphase (period 0.479 s vs 0.443 uncoupled).
#
# Adds to BOTH asim + aproj:
#   "c1 Commissural Inhibit"  Equil -70, SynAmp 2.749 (same numbers as RG
#                             Inhibit; dedicated type for GUI clarity)
#   "V3 Commissural Excite"   Equil -40, SynAmp 0.1   (weak; prevents latch)
# and repoints links cafe0023/25 (c1 out) and cafe0027/29 (V3 out).
# Usage: perl patch_comm_types.pl <dir>
use strict; use warnings;
my $dir = shift @ARGV or die "usage: patch_comm_types.pl <dir>\n";
my $ASIM = "$dir/Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim";
my $PROJ = "$dir/Walker_2_Layer_CPG_BilateralRG.aproj";
sub slurp { open my $f,'<',$_[0] or die; local $/; my $s=<$f>; close $f; $s }
sub spew { open my $o,'>',$_[0] or die; print {$o} $_[1]; close $o }
sub fg { my $n=shift; sprintf("cafe%04d-0000-4000-8000-%012d",$n,900000000000+$n) }

my $T_C1 = "cafe0041-0000-4000-8000-900000000041";
my $T_V3 = "cafe0042-0000-4000-8000-900000000042";
my $RGINH = "30b48a1d-6960-4e16-b812-a014fab17370";   # RG Inhibit
my $RGEXC = "9baf7f2a-bac2-4c80-b398-9dce2ab81f03";   # RG Excite
my $R2PF  = "bdb7ece4-859e-47cb-9d47-2e8156015e08";   # RG to PF Excite (SynAmp 0.1)

# ---- ASIM ----
{
  my $x = slurp($ASIM);
  if (index($x,"V3 Commissural Excite")>=0) { print "asim: already patched, skip\n"; }
  else {
  # new types after the RG to PF Excite type block
  my $tp = index($x,"<Name>RG to PF Excite</Name>");
  die "asim: RG to PF Excite type not found\n" if $tp<0;
  my $ts = rindex($x,"<SynapseType>",$tp); my $te = index($x,"</SynapseType>",$tp)+15;
  my $tmpl = substr($x,$ts,$te-$ts);
  my $c1type = $tmpl;
  $c1type =~ s/<Name>[^<]*<\/Name>/<Name>c1 Commissural Inhibit<\/Name>/;
  $c1type =~ s/<ID>[^<]*<\/ID>/<ID>$T_C1<\/ID>/;
  $c1type =~ s/<Equil>[^<]*<\/Equil>/<Equil>-70<\/Equil>/;
  $c1type =~ s/<SynAmp>[^<]*<\/SynAmp>/<SynAmp>2.749<\/SynAmp>/;
  my $v3type = $tmpl;
  $v3type =~ s/<Name>[^<]*<\/Name>/<Name>V3 Commissural Excite<\/Name>/;
  $v3type =~ s/<ID>[^<]*<\/ID>/<ID>$T_V3<\/ID>/;
  $v3type =~ s/<Equil>[^<]*<\/Equil>/<Equil>-40<\/Equil>/;
  $v3type =~ s/<SynAmp>[^<]*<\/SynAmp>/<SynAmp>0.1<\/SynAmp>/;
  substr($x,$te,0) = "\n".$c1type."\n".$v3type;
  # repoint links
  for my $n (23,25,27,29) {
    my $p = index($x,"<ID>".fg($n)."</ID>"); die "asim: link $n missing\n" if $p<0;
    my $s = rindex($x,"<Connexion>",$p); my $e = index($x,"</Connexion>",$p)+12;
    my $blk = substr($x,$s,$e-$s); my $nb = $blk;
    my $newtype = ($n==23||$n==25) ? $T_C1 : $T_V3;
    $nb =~ s{<SynapseTypeID>[^<]*<\/SynapseTypeID>}{<SynapseTypeID>$newtype</SynapseTypeID>} or die "asim: typeid swap $n\n";
    substr($x,$s,$e-$s) = $nb;
    print "asim: link $n -> ",(($n==23||$n==25)?"c1":"V3")," type\n";
  }
  spew($ASIM,$x); print "asim written\n";
  }
}

# ---- APROJ ----
{
  my $x = slurp($PROJ);
  if (index($x,"V3 Commissural Excite")>=0) { print "aproj: already patched, skip\n"; }
  else {
  # find the aproj synapse-type definition container: the Link-based types.
  # Aproj types are <Link> blocks with SynapseTypes.* ClassName and empty
  # Origin/Destination. Locate via a Link block whose Name is RG to PF Excite.
  my ($tblk);
  pos($x)=0;
  while ($x =~ /<Link>(.*?)<\/Link>/gs) {
    my $b=$1; next unless $b =~ /SynapseTypes\./;
    my ($nm)=$b=~/<Name>([^<]*)<\/Name>/;
    if (defined $nm && $nm eq 'RG to PF Excite') { $tblk="<Link>$b<\/Link>"; last; }
  }
  die "aproj: RG to PF Excite type Link not found\n" unless $tblk;
  my $c1type = $tblk;
  $c1type =~ s/<Name>[^<]*<\/Name>/<Name>c1 Commissural Inhibit<\/Name>/;
  $c1type =~ s/<ID>[^<]*<\/ID>/<ID>$T_C1<\/ID>/;
  $c1type =~ s{<EquilibriumPotential Value="[^"]*" Scale="milli" Actual="[^"]*"/>}{<EquilibriumPotential Value="-70" Scale="milli" Actual="-0.07"/>} or die "equil patch c1\n";
  $c1type =~ s{<MaxSynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>}{<MaxSynapticConductance Value="2.749" Scale="micro" Actual="2.749e-007"/>} or die "cond patch c1\n";
  my $v3type = $tblk;
  $v3type =~ s/<Name>[^<]*<\/Name>/<Name>V3 Commissural Excite<\/Name>/;
  $v3type =~ s/<ID>[^<]*<\/ID>/<ID>$T_V3<\/ID>/;
  $v3type =~ s{<EquilibriumPotential Value="[^"]*" Scale="milli" Actual="[^"]*"/>}{<EquilibriumPotential Value="-40" Scale="milli" Actual="-0.04"/>} or die "equil patch v3\n";
  $v3type =~ s{<MaxSynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>}{<MaxSynapticConductance Value="0.1" Scale="micro" Actual="1e-007"/>} or die "cond patch v3\n";
  substr($x,index($x,$tblk)+length($tblk),0) = "\n".$c1type."\n".$v3type;
  # repoint the 4 functional links + set their conductance
  for my $n (23,25,27,29) {
    my $p = index($x,"<ID>".fg($n)."</ID>"); die "aproj: link $n missing\n" if $p<0;
    my $s = rindex($x,"<Link>",$p); my $e = index($x,"</Link>",$p)+7;
    my $blk = substr($x,$s,$e-$s); my $nb = $blk;
    die "aproj: link $n not a synapse\n" unless $nb =~ /Behavior\.Synapse</;
    my $newtype = ($n==23||$n==25) ? $T_C1 : $T_V3;
    $nb =~ s{<SynapticTypeID>[^<]*<\/SynapticTypeID>}{<SynapticTypeID>$newtype</SynapticTypeID>} or die "aproj: typeid swap $n\n";
    my $cond = ($n==23||$n==25) ? "2.749" : "0.1";
    my $act  = ($n==23||$n==25) ? "2.749e-007" : "1e-008";
    $nb =~ s{<SynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>}{<SynapticConductance Value="$cond" Scale="micro" Actual="$act"/>} or die "aproj: cond swap $n\n";
    substr($x,$s,$e-$s) = $nb;
    print "aproj: link $n -> ",(($n==23||$n==25)?"c1":"V3")," type\n";
  }
  spew($PROJ,$x); print "aproj written\n";
  }
}
print "PATCH COMM TYPES DONE\n";
