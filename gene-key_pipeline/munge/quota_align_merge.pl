#!/usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
use Getopt::Long;

#python2 quota-alignment/cluster_utils.py --format=dag --log_evalue last/b73-phb47/b73-phb47.filtered.dag.go.aligncoords last/b73-phb47/b73-phb47.filtered.dag.go.aligncoords.qa

our ($infile, $outfile, $max_distance, $P, $QUOTA_ALIGN, $CLUSTER_UTILS,
    $CONFIG);

GetOptions(
    "infile|i=s"       => \$infile,
    "outfile|o=s"      => \$outfile,
    "max_distance|d=s" => \$max_distance);

run_quota_align_merge(
    infile   => $infile,
    outfile  => $outfile,
    max_dist => $max_distance);

sub run_quota_align_merge
{
    my %opts     = @_;
    my $infile   = $opts{infile};
    my $max_dist = $opts{max_dist};
    my $outfile  = $opts{outfile};

    #convert to quota-align format

    my $cmd = "quota-alignment/cluster_utils.py --format=dag --log_evalue $infile $infile.Dm$max_dist.qa";
    say "Converting dag output to quota_align format: $cmd";
    `$cmd`;

    $cmd = "quota-alignment/quota_align.py" . " --Dm=$max_dist --merge $infile.Dm$max_dist.qa";
    say "Running quota_align to merge diagonals:\n\t$cmd";
    `$cmd`;

    if (-r "$infile.Dm$max_dist.qa.merged") {
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
        open(IN,  "$infile.Dm$max_dist.qa.merged");
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
        say "The merged file $infile.Dm$max_dist.qa.merged was not created.";
    }
}
