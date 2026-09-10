#!/usr/bin/perl
# Compare element SHAPES (tagname + sorted attribute names) between pristine and current pages.
use strict; use warnings;
my ($origF, $curF) = @ARGV;

sub load_pages {
    my $f = shift;
    open my $fh, '<', $f or die; local $/; my $x = <$fh>; close $fh;
    my %pages;
    pos($x) = 0;
    while ($x =~ /<DiagramXml><!\[CDATA\[(.*?)\]\]><\/DiagramXml>/gs) {
        my $d = $1;
        my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
        $pn =~ s/&amp;/&/g;
        $pages{$pn} = $d;
    }
    return %pages;
}
sub shapes {
    my $d = shift;
    my %shapes;   # "tag a,b,c" -> count
    while ($d =~ /<([A-Za-z]+)((?:\s+[-\w]+="[^"]*")*)\s*(\/?)>/g) {
        my ($tag, $attrs, $selfclose) = ($1, $2, $3);
        my @names = $attrs =~ /([-\w]+)=/g;
        $shapes{"$tag [" . join(',', sort @names) . "]"}++;
    }
    return \%shapes;
}

my %o = load_pages($origF);
my %c = load_pages($curF);
for my $pn (sort keys %c) {
    my $os = $o{$pn} ? shapes($o{$pn}) : {};
    my $cs = shapes($c{$pn});
    my @new = grep { !$os->{$_} } sort keys %$cs;
    my @gone = grep { !$cs->{$_} } sort keys %$os;
    next unless @new || @gone;
    print "== $pn ==\n";
    for (@new)  { print "  + $_  (x$cs->{$_})\n"; }
    for (@gone) { print "  - $_  (x$os->{$_})\n"; }
}
print "done\n";
