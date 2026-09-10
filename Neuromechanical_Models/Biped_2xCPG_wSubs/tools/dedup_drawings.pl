#!/usr/bin/perl
# Remove duplicate drawing links: on each page, if a b0a0/b1a0 drawing link Tag appears
# more than once, keep the LAST occurrence (the corrected redraw) and drop earlier ones.
use strict; use warnings;
my ($in, $out) = @ARGV;
open my $fh, '<', $in or die $!; local $/; my $xml = <$fh>; close $fh;

my $result = '';
my $pos = 0;
my $removed = 0;
while ($xml =~ /<DiagramXml><!\[CDATA\[/g) {
    my $cdata_start = pos($xml);
    my $cdata_end = index($xml, ']]></DiagramXml>', $cdata_start);
    last if $cdata_end < 0;
    my $d = substr($xml, $cdata_start, $cdata_end - $cdata_start);
    my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
    $pn =~ s/&amp;/&/g;

    my $af = index($d, '<AddFlow ');
    my $body_start = index($d, '>', $af) + 1;
    my $epi = rindex($d, '</AddFlow>');
    my $body = substr($d, $body_start, $epi - $body_start);

    # count occurrences of each drawing-link Tag to know which are duplicates
    my %count;
    pos($body) = 0;
    while ($body =~ /<Link\b[^>]*>(.*?)<\/Link>/gs) {
        my ($t) = $1 =~ /<Tag>([^<]+)</;
        $count{$t}++ if $t && $t =~ /^b[01]a0/;
    }

    # rebuild body: skip earlier duplicate occurrences (keep last), skip 1st+? NO: keep last
    my %seen;
    my $newbody = '';
    my $last = 0;
    my $removedHere = 0;
    pos($body) = 0;
    while ($body =~ /(<Link\b[^>]*>)(.*?)(<\/Link>)/gs) {
        my ($open, $inner, $close) = ($1, $2, $3);
        my $mstart = $-[0];
        my $mend = $+[0];
        my ($t) = $inner =~ /<Tag>([^<]+)</;
        my $isDup = (defined $t && $t =~ /^b[01]a0/ && $count{$t} > 1);
        $newbody .= substr($body, $last, $mstart - $last);
        if ($isDup) {
            $seen{$t}++;
            if ($seen{$t} < $count{$t}) {
                # earlier duplicate: drop
                $removedHere++;
                $last = $mend;
                next;
            }
        }
        $newbody .= substr($body, $mstart, $mend - $mstart);
        $last = $mend;
    }
    $newbody .= substr($body, $last);
    if ($removedHere) {
        $newbody =~ s{(<AddFlow Nodes="(\d+)" Links=")(\d+)(")}{$1 . ($3 - $removedHere) . $4}e;
        substr($d, $body_start, $epi - $body_start) = $newbody;
    }
    $removed += $removedHere;
    print "page '$pn': removed $removedHere duplicate drawings\n" if $removedHere;
    $result .= substr($xml, $pos, $cdata_start - $pos);
    $result .= $d;
    $pos = $cdata_end;
}
$result .= substr($xml, $pos);
print "SUMMARY: removed $removed duplicate drawings\n";
open my $ofh, '>', $out or die $!;
print $ofh $result;
close $ofh;
print "WROTE $out\n";
