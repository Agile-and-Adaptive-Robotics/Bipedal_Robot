#!/usr/bin/perl
# sweep_comm.pl — air-gain sweep for commissural pathways on the ASIM only.
# Patches <G> of c1 outputs (cafe0023/25) and V3 outputs (cafe0027/29),
# runs AnimatSimulator in a scratch dir, scores oscillation.
# Usage: perl sweep_comm.pl <asim> <outdir>
use strict; use warnings;
my ($asim,$outdir) = @ARGV;
my $SIMEXE = "D:/Program Files (x86)/NeuroRobotic Technologies/AnimatLab/bin/AnimatSimulator.exe";

sub slurp { open my $f,'<',$_[0] or die; local $/; my $s=<$f>; close $f; $s }
sub patchG {
  my ($xr,$idsuffix,$g)=@_;
  my $tag = "cafe$idsuffix";
  my $pos=0;
  while (($pos=index($$xr,"<ID>$tag",$pos))>=0) {
    my $s=rindex($$xr,'<Connexion>',$pos); my $e=index($$xr,'</Connexion>',$pos)+12;
    my $blk=substr($$xr,$s,$e-$s);
    if ($blk =~ /^<Connexion>\s*<ID>cafe\d{4}-/) {
      my $nb=$blk; $nb =~ s/<G>[^<]*<\/G>/<G>$g<\/G>/;
      substr($$xr,$s,$e-$s)=$nb; return 1;
    }
    $pos+=8;
  }
  return 0;
}

sub score {
  my ($txt)=@_;
  open my $f,'<',$txt or return {err=>"no chart"};
  my $h=<$f>; $h=~s/\r//g; my @n=split/\t/,$h;
  my %c; for my $i (0..$#n){ $c{$n[$i]}=$i; }
  my %on; my %mn; my %mx; my @prev;
  while(<$f>){ s/\r//g; my @v=split/\t/;
    for my $nm ("L RG ext","R RG ext") {
      my $i=$c{$nm};
      $mn{$nm}=$v[$i] if !defined($mn{$nm})||$v[$i]<$mn{$nm};
      $mx{$nm}=$v[$i] if !defined($mx{$nm})||$v[$i]>$mx{$nm};
      if (defined $prev[$i]) { push @{$on{$nm}},$v[1] if $prev[$i]<=-0.060 && $v[$i]>-0.060; }
      $prev[$i]=$v[$i];
    }
  }
  close $f;
  my %r;
  for my $nm ("L RG ext","R RG ext") {
    my $o=$on{$nm}//[];
    $r{"n$nm"}=scalar @$o;
    if (@$o>=3) { my @p=map{$o->[$_]-$o->[$_ - 1]} 1..$#$o; my $s=0; $s+=$_ for @p; $r{"T$nm"}=sprintf("%.3f",$s/@p); }
  }
  # antiphase: interleave of onsets
  if (($r{"nL RG ext"}//0)>=3 && ($r{"nR RG ext"}//0)>=3) {
    my @all = sort { $a<=>$b } (@{$on{"L RG ext"}}, @{$on{"R RG ext"}});
    my $inter=1;
    for my $i (0..$#all-1) { $inter=0 if abs($all[$i+1]-$all[$i])<0.02; }
    $r{interleave}=$inter;
  }
  $r{lock} = (($r{"nL RG ext"}//0)<2 && ($r{"nR RG ext"}//0)<2) ? 1 : 0;
  return \%r;
}

my @c1g = (0.05, 0.1, 0.2, 0.5);
my @v3g = (0.01, 0.05, 0.15, 0.5);
mkdir $outdir unless -d $outdir;
printf("%-6s %-6s | %-14s %-14s %s\n","c1out","v3out","L ext","R ext","interleave");
for my $c (@c1g) {
  for my $v (@v3g) {
    my $x = slurp($asim);
    patchG(\$x,"0023",$c) or die "c1 link 23 missing";
    patchG(\$x,"0025",$c) or die "c1 link 25 missing";
    patchG(\$x,"0027",$v) or die "v3 link 27 missing";
    patchG(\$x,"0029",$v) or die "v3 link 29 missing";
    my $dir = "$outdir/c$c" . "_v$v"; mkdir $dir;
    my $tmp = "$dir/test.asim";
    open my $o,'>',$tmp; print $o $x; close $o;
    system("\"$SIMEXE\" \"$tmp\" > $dir/run.log 2>&1");
    my $r = score("$dir/Rhythm Generator.txt");
    my $ls = sprintf("%db%s",$r->{"nL RG ext"}//0, $r->{"T L RG ext"} ? " T=".$r->{"T L RG ext"} : "");
    my $rs = sprintf("%db%s",$r->{"nR RG ext"}//0, $r->{"T R RG ext"} ? " T=".$r->{"T R RG ext"} : "");
    printf("%-6s %-6s | %-14s %-14s %s%s\n",$c,$v,$ls,$rs, ($r->{interleave}?"YES":""), ($r->{lock}?" LOCK":""));
  }
}
print "sweep done\n";
