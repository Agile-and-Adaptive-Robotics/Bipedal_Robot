#!/usr/bin/perl
# Set each CDATA page's declared AddFlow Nodes/Links counts to the actual element counts.
use strict; use warnings;
my ($in, $out) = @ARGV;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

my $result = '';
my $pos = 0;
my $fixed = 0;
while ($xml =~ /<DiagramXml><!\[CDATA\[/g) {
    my $s = pos($xml);
    my $e = index($xml, ']]></DiagramXml>', $s);
    last if $e < 0;
    my $d = substr($xml, $s, $e - $s);
    my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
    $pn =~ s/&amp;/&/g;

    my ($na, $la) = $d =~ /<AddFlow Nodes="(\d+)" Links="(\d+)"/;
    my $af = index($d, '<AddFlow ');
    my $body_start = index($d, '>', $af) + 1;
    my $epi = rindex($d, '</AddFlow>');
    my $body = substr($d, $body_start, $epi - $body_start);
    my ($no, $lo) = (0, 0);
    pos($body) = 0;
    while ($body =~ /<(Node|Link)\b/g) {
        my $t = $1;
        my $cs = pos($body);
        my $closeTag = $t eq 'Node' ? '</Node>' : '</Link>';
        my $ce = index($body, $closeTag, $cs);
        last if $ce < 0;
        $no++ if $t eq 'Node'; $lo++ if $t eq 'Link';
        pos($body) = $ce + length($closeTag);
    }
    if ($na != $no || $la != $lo) {
        my ($openTag) = $d =~ /(<AddFlow \b[^>]*>)/;
        my $newOpen = $openTag;
        $newOpen =~ s/Nodes="\d+"/Nodes="$no"/;
        $newOpen =~ s/Links="\d+"/Links="$lo"/;
        my $opos = index($d, $openTag);
        substr($d, $opos, length($openTag)) = $newOpen;
        substr($xml, $s, $e - $s) = $d;
        $fixed++;
        print "$pn: declared N=$na L=$la -> N=$no L=$lo\n";
    }
    pos($xml) = $s;
}
open my $ofh, '>', $out or die $!;
print $ofh $xml;
close $ofh;
print "fixed $fixed pages\n";
