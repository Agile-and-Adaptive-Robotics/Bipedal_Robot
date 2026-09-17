#!/usr/bin/perl
# verify_handles.pl <aproj> — test how AddFlow drawing Link Org/Dst encode endpoints.
# Hypothesis: Org/Dst = position (0- or 1-based) of the endpoint NODE ENTRY in the
# INTERLEAVED file order (nodes+links counted together). Alternatives: node-only
# ordinal. Tests against every original drawing whose Tag matches a functional link.
use strict; use warnings;
my $f = shift @ARGV;
local $/; open my $fh,'<',$f or die; my $x=<$fh>; close $fh;
my ($cd)=$x=~/<DiagramXml><!\[CDATA\[(.*?)\]\]><\/DiagramXml>/s or die "no cdata";

# ordered shape list
my @shapes;   # [kind, tag]
my $pos=0;
while ($cd =~ /(<Node\s+Left=[^>]*>|<Link\s+Org=[^>]*>)/g) {
  my $open=$1;
  my $kind = $open=~/^<Node/ ? "N" : "L";
  my $from = pos($cd);
  my $close = $kind eq "N" ? "</Node>" : "</Link>";
  my $e = index($cd,$close,$from);
  my $blk = substr($cd,$from,$e-$from);
  my ($tag) = $blk =~ /<Tag>([^<]*)<\/Tag>/;
  push @shapes,[$kind,$tag//""];

}
print "shapes: ",scalar(@shapes)," (N=",scalar(grep{$_->[0] eq 'N'}@shapes)," L=",scalar(grep{$_->[0] eq 'L'}@shapes),")\n";
my %idx0; my %idxN; my $n=0;
for my $i (0..$#shapes) { $idx0{$shapes[$i][1]}=$i if $shapes[$i][1] && !exists $idx0{$shapes[$i][1]};
  if ($shapes[$i][0] eq 'N') { $idxN{$shapes[$i][1]}=$n++ if $shapes[$i][1] && !exists $idxN{$shapes[$i][1]}; } }

# functional links (Behavior.Synapse + Behavior.Links.Adapter)
my %func;   # linkid => [origin,dest,isAdapter]
pos($x)=0;
while ($x =~ /<Link>(.*?)<\/Link>/gs) {
  my $b=$1;
  my $isAd = $b=~/Behavior\.Links\.Adapter</ ? 1 : 0;
  next unless $isAd || $b=~/Behavior\.Synapse</;
  my ($id)=$b=~/<ID>([^<]+)<\/ID>/; my ($o)=$b=~/<OriginID>([^<]+)<\/OriginID>/; my ($d)=$b=~/<DestinationID>([^<]+)<\/DestinationID>/;
  next unless $id && $o && $d;
  $func{$id}=[$o,$d,$isAd];
}
# which functional links are DRAWN?
my $drawn=0; my $undrawn=0; my @undrawnAd; my @undrawnSyn;
for my $id (keys %func) {
  if (exists $idx0{$id}) { $drawn++; }
  else { $undrawn++; push @{$func{$id}[2]?\@undrawnAd:\@undrawnSyn},$id; }
}
print "functional links: ",scalar(keys %func),"  drawn=$drawn undrawn=$undrawn (undrawn adapters=",scalar(@undrawnAd)," synapses=",scalar(@undrawnSyn),")\n";

# convention test on drawn SYNAPSE links (neural endpoints, most reliable)
my (%hit0,%hit1,%hitN,%tot);
pos($cd)=0;
while ($cd =~ /<Link\s+Org="(\d+)"\s+Dst="(\d+)">/g) {
  my ($org,$dst)=($1,$2);
  my $from=pos($cd); my $e=index($cd,"</Link>",$from);
  my ($tag)=substr($cd,$from,$e-$from)=~/<Tag>([^<]*)<\/Tag>/;
  next unless $tag && $func{$tag};
  my ($o,$d,$ad)=@{$func{$tag}};
  next if $ad;   # adapter links may involve undrawn physicals; test separately
  my $o0=$idx0{$o}//-1; my $d0=$idx0{$d}//-1;
  my $oN=$idxN{$o}//-1; my $dN=$idxN{$d}//-1;
  $tot{all}++;
  $hit0{all}++ if $org==$o0 && $dst==$d0;
  $hit1{all}++ if $org==$o0+1 && $dst==$d0+1;
  $hitN{all}++ if ($org==$oN||$org==$oN+1) && ($dst==$dN||$dst==$dN+1);
  # partial: org only
  $hit0{org}++ if $org==$o0; $hit1{org}++ if $org==$o0+1;
  $tot{org}++;
}
print "interleaved 0-based: full=$hit0{all}/$tot{all}  orgOnly=$hit0{org}/$tot{org}\n";
print "interleaved 1-based: full=$hit1{all}/$tot{all}  orgOnly=$hit1{org}/$tot{org}\n";
print "node-ordinal(0/1):   full=$hitN{all}/$tot{all}\n";

# adapter links: how many drawn, and do their endpoints resolve?
my ($adDr,$adOk0)=(0,0);
pos($cd)=0;
while ($cd =~ /<Link\s+Org="(\d+)"\s+Dst="(\d+)">/g) {
  my ($org,$dst)=($1,$2);
  my $from=pos($cd); my $e=index($cd,"</Link>",$from);
  my ($tag)=substr($cd,$from,$e-$from)=~/<Tag>([^<]*)<\/Tag>/;
  next unless $tag && $func{$tag} && $func{$tag}[2];
  $adDr++;
  my ($o,$d)=@{$func{$tag}}[0,1];
  $adOk0++ if ($idx0{$o}//-2)==$org && ($idx0{$d}//-2)==$dst;
}
print "adapter drawings: $adDr (endpoint-exact 0-based: $adOk0)\n";
# do physical bodies appear as drawn node tags? sample: first adapter's OriginID
my $sample = (grep { $_->[2] } values %func)[0];
if ($sample) {
  print "sample adapter link origin=$sample->[0] drawnAsNode=",($idx0{$sample->[0]}//":")," dest=$sample->[1] drawnAsNode=",($idx0{$sample->[1]}//":"),"\n";
}
