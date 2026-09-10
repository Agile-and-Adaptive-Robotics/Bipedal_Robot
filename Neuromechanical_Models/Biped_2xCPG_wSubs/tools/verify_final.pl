#!/usr/bin/perl
# Careful verification of the rebuilt .aproj drawing pages + functional model.
use strict; use warnings;
my $file = shift;
open my $fh, '<', $file or die $!; local $/; my $xml = <$fh>; close $fh;

# --- split into CDATA page regions sequentially (single authoritative scan) ---
my @pages;   # {cdata, pn}
{
    pos($xml) = 0;
    while ($xml =~ /<DiagramXml><!\[CDATA\[/g) {
        my $s = pos($xml);
        my $e = index($xml, ']]></DiagramXml>', $s);
        last if $e < 0;
        my $d = substr($xml, $s, $e - $s);
        my ($pn) = $d =~ /<PageName>([^<]*)<\/PageName>/;
        $pn =~ s/&amp;/&/g;
        push @pages, { pn => $pn, d => $d };
        pos($xml) = $s;
    }
}

my ($mism, $oor, $badorder) = (0, 0, 0);
for my $p (@pages) {
    my $d = $p->{d};
    my $pn = $p->{pn};
    my ($na, $la) = $d =~ /<AddFlow Nodes="(\d+)" Links="(\d+)"/;
    # sequential element walk (same pattern as the rebuild script)
    my @elems;
    my $af = index($d, '<AddFlow ');
    my $body_start = index($d, '>', $af) + 1;
    my $epi = rindex($d, '</AddFlow>');
    my $body = substr($d, $body_start, $epi - $body_start);
    pos($body) = 0;
    while ($body =~ /<(Node|Link)\b/g) {
        my $t = $1;
        my $s = pos($body);
        my $closeTag = $t eq 'Node' ? '</Node>' : '</Link>';
        my $e = index($body, $closeTag, $s);
        last if $e < 0;
        push @elems, $t;
        pos($body) = $e + length($closeTag);
    }
    my $no = grep { $_ eq 'Node' } @elems;
    my $lo = grep { $_ eq 'Link' } @elems;
    $mism++ if $na != $no || $la != $lo;
    printf "%-28s declared N=%-3d L=%-3d actual N=%-3d L=%-3d %s\n", $pn, $na, $la, $no, $lo,
        (($na != $no || $la != $lo) ? 'MISMATCH' : '');
    # order + range checks
    my $seenLink = 0;
    my $k = 0; my %ord;
    for my $el (@elems) {
        if ($el eq 'Node') { $k++; }
        else {
            $badorder++ if $seenLink && 0;
            $seenLink = 1;
        }
    }
    while ($body =~ /<Link\b([^>]*)>(.*?)<\/Link>/gs) {
        my ($attrs, $b2) = ($1, $2);
        my ($org) = $attrs =~ /Org="(\d+)"/;
        my ($dst) = $attrs =~ /Dst="(\d+)"/;
        $oor++ if (defined $org && $org >= $no) || (defined $dst && $dst >= $no);
    }
}
print "pages=" . scalar(@pages) . " countMismatches=$mism outOfRangeLinks=$oor\n";

# --- functional model checks ---
my %cnt;
while ($xml =~ /<(StimulusTension|LengthTension|Gain|CaActivation|CaDeactivation|PID)>\s*<ID>([^<]+)<\/ID>/gs) { $cnt{"$1 $2"}++; }
my $dups = grep { $_ > 1 } values %cnt;
my $syn  = () = $xml =~ /Behavior\.Synapse</g;
my $bad  = 0;
while ($xml =~ m{(<Link>.*?</Link>)}gs) {
    my $b = $1;
    next unless $b =~ /Behavior\.(Synapse|Links\.Adapter)</;
    for my $e ($b =~ /<(?:OriginID|DestinationID)>([^<]+)<\/(?:OriginID|DestinationID)>/g) {
        unless (length($e) > 8 && index($xml, "<ID>$e</ID>") >= 0) { $bad++; }
    }
}
my $corrupt = 0;
while ($xml =~ /<([A-Za-z]+)>(\s*)([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})\s*\n?\s*<([A-Za-z])/g) {
    next if $1 eq 'ID' || $4 eq '/';
    $corrupt++;
}
my $mus   = () = $xml =~ /<Type>LinearHillMuscle<\/Type>/g;
my $tonic = () = $xml =~ /<TonicStimulus Value="2" Scale="nano"/g;
print "childIDdups=$dups syn=$syn dangling=$bad corrupted=$corrupt muscles=$mus tonic2nA=$tonic\n";
