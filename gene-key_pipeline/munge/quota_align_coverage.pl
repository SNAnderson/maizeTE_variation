#!/usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
use Getopt::Long;

our (
    $infile,      $outfile,       $orgratio1,
    $orgratio2,   $overlap,       $prefix);

GetOptions(
    "infile|if=s"           => \$infile,
    "outfile|of=s"          => \$outfile,
    "depth_ratio_org1|d1=s" => \$orgratio1,
    "depth_ratio_org2|d2=s" => \$orgratio2,
    "depth_overlap|o=s"     => \$overlap,
    "prefix|p=s"            => \$prefix);


run_quota_align_coverage(
    infile  => $infile,
    outfile => $outfile,
    org1    => $orgratio1,
    org2    => $orgratio2,
    overlap => $overlap,
    prefix  => $prefix);

sub run_quota_align_coverage
{
    my %opts         = @_;
    my $infile       = $opts{infile};
    my $org1         = $opts{org1};      #ratio of org1
    my $org2         = $opts{org2};      #ratio of org2
    my $overlap_dist = $opts{overlap};
    my $outfile      = $opts{outfile};
    my $prefix       = $opts{prefix};

    #convert to quota-align format
    my $cov_cmd = "quota-alignment/cluster_utils.py --format=dag --log_evalue $infile $prefix";
    my $qa_cmd = "quota-alignment/quota_align.py --Nm=$overlap_dist --quota=$org1:$org2 $prefix";

    say "Convert command: $cov_cmd";
    say "Quota Align command: $qa_cmd";

    say "Converting dag output to quota_align format.";
    `$cov_cmd`;
    say "Running quota_align to find syntenic coverage.";
    `$qa_cmd`;

    if (-r "$prefix") {
        my %data;
        $/ = "\n";
        open(IN, $infile);
        while (<IN>) {
            next if /^#/;
            my @line = split /\t/;
            $data{join("_", $line[0], $line[2], $line[4], $line[6])} = $_;
        }
        close IN;
        open(OUT, ">$outfile");
        open(IN,  "$prefix");
        while (<IN>) {
            if (/^#/) {
                print OUT $_;
            } else {
                chomp;
                my @line = split /\t/;
                print OUT $data{
                    join("_", $line[0], $line[1], $line[2], $line[3])};
            }
        }
        close IN;
        close OUT;
    } else {
        say "Syntenic coverage failed to output $prefix";
    }
}
