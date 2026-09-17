#!/usr/bin/perl
# ground_test.pl — build a GROUND-TEST asim variant in a scratch dir.
#   - Root: Freeze -> False, Position.y -> $rooty
#   - adds a "Contact" FileChart: ContactCount of toe_L/toe_R/foot_L/foot_R
#     contact bodies (diagnostic for stance/tuning)
# Usage: perl ground_test.pl <src.asim> <outdir> <rooty> [shelfy]
#   shelfy: if given, also adds a static pelvis-shelf box (partial support)
#           at y=shelfy under the Root (x/z from Root), size 0.30x0.04x0.40 m.
use strict; use warnings;
my ($src,$outdir,$rooty,$shelfy) = @ARGV;
die "usage: ground_test.pl <src.asim> <outdir> <rooty> [shelfy]\n" unless $src && $outdir && defined $rooty;

sub slurp { open my $f,'<',$src or die; local $/; my $s=<$f>; close $f; $s }
my $x = slurp($src);

# contact body ids by name
my %cid;
for my $bn ("toe_L_contact","toe_R_contact","foot_L_contact","foot_R_contact") {
  my $i=index($x,"<Name>$bn</Name>");
  die "body $bn not found\n" if $i<0;
  my $s=rindex($x,"<RigidBody>",$i);
  my ($id) = substr($x,$i,300) =~ /<ID>([^<]+)<\/ID>/;
  die "id for $bn\n" unless $id;
  $cid{$bn}=$id;
}

# --- Root: unfreeze + lower ---
{
  my $i=index($x,"<Name>Root</Name>");
  die "Root not found\n" if $i<0;
  my $s=rindex($x,"<RigidBody>",$i); my $e=index($x,"</RigidBody>",$i);
  my $blk=substr($x,$s,$e-$s);
  my $nb=$blk;
  $nb =~ s/<Freeze>True<\/Freeze>/<Freeze>False<\/Freeze>/ or die "Root freeze\n";
  $nb =~ s/(<Position [^>]*?)y="[-\d.e]+"([^>]*\/>)/${1}y="$rooty"$2/s or die "Root y\n";
  substr($x,$s,$e-$s)=$nb;
  print "Root: freeze=False y=$rooty\n";
}

# --- optional pelvis shelf (clone WalkingPath box, top-level static) ---
if (defined $shelfy) {
  my $i=index($x,"<Name>WalkingPath</Name>");
  die "WalkingPath not found\n" if $i<0;
  my $s=rindex($x,"<RigidBody>",$i); my $e=index($x,"</RigidBody>",$i)+12;
  my $blk=substr($x,$s,$e-$s);
  my $nb=$blk;
  my $gseq = 950;
  my $ng = sub { sprintf("cafe%04d-0000-4000-8000-%012d",$gseq,900000000000+$gseq++) };
  my $nid = $ng->();
  $nb =~ s/<ID>[^<]+<\/ID>/<ID>$nid<\/ID>/ or die;
  $nb =~ s/<Name>WalkingPath<\/Name>/<Name>PelvisShelf<\/Name>/ or die;
  $nb =~ s/(<Position [^>]*?)x="[-\d.e]+"([^>]*\/>)/${1}x="-3.454"$2/s or die;   # under Root
  $nb =~ s/(<Position [^>]*?)y="[-\d.e]+"([^>]*\/>)/${1}y="$shelfy"$2/s or die;
  $nb =~ s/(<Position [^>]*?)z="[-\d.e]+"([^>]*\/>)/${1}z="0"$2/s or die;
  $nb =~ s/<Length>[\d.e+-]+<\/Length>/<Length>0.30<\/Length>/ or die;
  $nb =~ s/<Width>[\d.e+-]+<\/Width>/<Width>0.40<\/Width>/ or die;
  $nb =~ s/<Height>[\d.e+-]+<\/Height>/<Height>0.04<\/Height>/ or die;
  $nb =~ s/<Mass>[\d.e+-]+<\/Mass>/<Mass>9e+008<\/Mass>/ or die;
  $nb =~ s/<Density>[\d.e+-]+<\/Density>/<Density>10<\/Density>/ or die;
  # insert after WalkingPath block (top-level sibling)
  substr($x,$e,0)="\n".$nb;
  print "PelvisShelf: top at y=",($shelfy+0.02),"\n";
}

# --- Contact chart: clone the Rhythm Generator chart, repoint ---
{
  my $i=index($x,"<OutputFilename>Rhythm Generator.txt<\/OutputFilename>");
  die "RG chart not found\n" if $i<0;
  my $cs=rindex($x,"<DataChart>",$i); my $ce=index($x,"</DataChart>",$i)+12;
  my $blk=substr($x,$cs,$ce-$cs);
  # strip its DataColumns
  my $nb = $blk;
  $nb =~ s/<DataColumns>.*?<\/DataColumns>/<DataColumns>PLACEHOLDER<\/DataColumns>/s or die "cols strip\n";
  my $gseq = 970;
  my $ng = sub { sprintf("cafe%04d-0000-4000-8000-%012d",$gseq,900000000000+$gseq++) };
  # column template from the (now removed) first column
  my ($colt) = $blk =~ /(<DataColumn>.*?<\/DataColumn>)/s;
  die "col template\n" unless $colt;
  my @names = ("toe_L_contact","toe_R_contact","foot_L_contact","foot_R_contact");
  my $cols='';
  for my $n (@names) {
    my $c=$colt;
    my ($oid)=$c=~/<ID>([^<]+)/;
    $c =~ s/<ID>\Q$oid\E<\/ID>/<ID>@{[$ng->()]}\n<\/ID>/ or die;
    $c =~ s/<ColumnName>[^<]*<\/ColumnName>/<ColumnName>$n<\/ColumnName>/ or die;
    $c =~ s/<DataType>[^<]*<\/DataType>/<DataType>ContactCount<\/DataType>/ or die;
    $c =~ s/<TargetID>[^<]*<\/TargetID>/<TargetID>$cid{$n}<\/TargetID>/ or die;
    $cols .= $c."\n";
  }
  $nb =~ s/PLACEHOLDER/$cols/;
  $nb =~ s/<OutputFilename>[^<]*<\/OutputFilename>/<OutputFilename>Contact.txt<\/OutputFilename>/ or die;
  $nb =~ s/<Name>[^<]*<\/Name>/<Name>Contact<\/Name>/ or die;
  my ($cid2)=$blk=~/<ID>([^<]+)/;
  $nb =~ s/<ID>\Q$cid2\E<\/ID>/<ID>@{[$ng->()]}\n<\/ID>/ or die;
  substr($x,$ce,0)="\n".$nb;
  print "Contact chart added (4 ContactCount columns)\n";
}

mkdir $outdir unless -d $outdir;
open my $o,'>',"$outdir/test.asim" or die; print {$o} $x; close $o;
print "wrote $outdir/test.asim\n";
