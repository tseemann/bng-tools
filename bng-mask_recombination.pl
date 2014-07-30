#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $symbol);
setOptions();

my($afa, $tab) = @ARGV;
-r $afa or die "Can not open alignment file: $afa";
-r $tab or die "Can not open bratNextGen tabular file: $tab";

#>BPH0530c
#gtgtcactttcgctttggcagcagtgtcttgcccgattgcaggatgagtt
#>EC_907700_uid183806
#gtgtcactttcgctttggcagcagtgtcttgcccgattgcaggatgagtt
#>EC_BIDMC_20A_uid202030
#gtgtcactttcgctttggcagcagtgtcttgcccgattgcaggatgagtt

my $bp=0;
my @id;
my %seq;
my $afa_fh = Bio::SeqIO->new(-file=>$afa, -format=>'fasta');
while (my $seq = $afa_fh->next_seq) {
  push @id, $seq->id;  # keep order
  print STDERR "Loading alignment $. $id[-1] ...";
  $seq{ $seq->id } = $seq->seq;  
  my $L = $seq->length;
  print STDERR " length $L\n";
  $bp += $L;
}

#LIST OF FOREIGN GENOMIC SEGMENTS:
#Start     End       Origin  HomeCluster  BAPSIndex  StrainName
#758790    759380    2       3            1          BPH0530c
#766054    766384    4       3            1          BPH0530c
#766054    766384    4       3            3          EC_907700_uid183806
# 85279    785990    2       3            10         EC958
#792281    794945    2       3            10         EC958

my $masked=0;

open my $tab_fh, '<', $tab;
while (<$tab_fh>) {
  chomp;
  my @x = split ' ';  # white space separated....
  next unless @x == 6;
  next unless $x[0] =~ m/^\d+$/;
  my $start = $x[0] - 1;
  my $len = $x[1] - $x[0] + 1;
  my $id = $x[5];
  substr $seq{$id}, $start, $len, ($symbol)x$len;
  $masked += $len;
}

# write out the new masked alignment
my $out_fh = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
for my $id (@id) {
  my $seq = Bio::Seq->new(-id=>$id, -seq=>$seq{$id});
  printf STDERR "Writing masked alignment: $id ... length %d\n", $seq->length;
  $out_fh->write_seq($seq );
}

printf STDERR "Masked %d of %d bases (%.2f%%)\n", $masked, $bp, ($masked*100/$bp);
print STDERR "Done.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"symbol=s",  VAR=>\$symbol, DEFAULT=>'?', DESC=>"Masking symbol"},
  );

  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] scapper.aln bratrecomb.tab > masked.aln\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
