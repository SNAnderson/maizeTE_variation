#!/usr/bin/perl -w
use strict;

open ( my $in, '<', $ARGV[0] ) or die;
open ( my $out, '>', $ARGV[1] ) or die;
open ( my $out2, '>', $ARGV[2] ) or die;

<$in>;
while ( my $line = <$in> ) {
	chomp $line;

	if ( $line !~ /^(\s+)?\d+/ ) {
		next;
	}

    my @fields;
	for my $field (split ' ', $line) {
        push @fields, $field if ( $field !~ /\|/ );
    }

	my $queryStart = $fields[0];
	my $queryStop = $fields[1];
	my $hitStart = $fields[2];
	my $hitStop = $fields[3];
	my $queryChr = $fields[11]; $queryChr =~ s/chr//;
	my $hitChr = $fields[12]; $hitChr =~ s/lcl\|//;

	print $out "$queryChr\t";
	$queryStart > $queryStop ? print $out $queryStop-1,"\t$queryStart\t" : print $out $queryStart-1,"\t$queryStop\t";
	$hitStart > $hitStop ? print $out "$hitChr-",$hitStop-1,";$hitStart" : print $out "$hitChr-",$hitStart-1,";$hitStop";
	print $out "\n";

  print $out2 "$hitChr\t";
	$hitStart > $hitStop ? print $out2 $hitStop-1,"\t$hitStart\t" : print $out2 $hitStart-1,"\t$hitStop\t";
	$queryStart > $queryStop ? print $out2 "$queryChr-",$queryStop-1,";$queryStart" : print $out2 "$queryChr-",$queryStart-1,";$queryStop";
	print $out2 "\n";
}
close $in;
close $out;
close $out2;
