#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
#use Bio::TreeIO;
#use Bio::AlignIO;
use Data::Dumper;
use List::Util qw(min max);
use List::MoreUtils qw(uniq);
use SVG;

my(@Options, $verbose, $tab_fn, $tree_fn, $aln_fn, $svg_fn,
             $width, $height, $consensus, $border, $monochrome);
setOptions();

my $aln_fh = get_fh($aln_fn, 'FASTA alignment');
my $tab_fh = get_fh($tab_fn, 'BRATNextGen tabular');
$svg_fn or die "Please provide an output filename eg. --svg out.svg";

#my $tree_fh = get_fh($tree_fn, 'Newick tree');

#my $aln_io = Bio::AlignIO->new(-fh=>$aln_fh, -format=>'fasta');
#my $aln = $aln_io->next_aln;
my $aln_io = Bio::SeqIO->new(-fh=>$aln_fh, -format=>'fasta');
my $aln = $aln_io->next_seq;
my $length = $aln->length;
print STDERR "Length of alignment: $length\n";

#LIST OF FOREIGN GENOMIC SEGMENTS:
#Start     End       Origin  HomeCluster  BAPSIndex  StrainName
#758790    759380    2       3            1          BPH0530c
#766054    766384    4       3            1          BPH0530c
#766054    766384    4       3            3          EC_907700_uid183806
# 85279    785990    2       3            10         EC958
#792281    794945    2       3            10         EC958

my %box;
my $boxes;
my %colour_of;
my %group_of;

while (<$tab_fh>) {
  chomp;
  my @x = split ' ';  # white space separated....
  next unless @x == 6;
  next unless $x[0] =~ m/^\d+$/;
  my $start = $x[0] - 1;
  my $len = $x[1] - $x[0] + 1;
  my $id = $x[5];
  my $group = $x[2];
  $colour_of{$group}++;
  my $home = $x[3];
  $group_of{$id} = $home;
  push @{$box{$id}}, [ $start, $len, $group ];
  $boxes++;
}
print STDERR "Need to draw $boxes boxes.\n";

print STDERR "Choosing colours...\n";
for my $n (uniq sort values %group_of) {
  my $colour = <DATA>;
  chomp $colour;
  $colour_of{$n} = $monochrome ? 'grey' : $colour;
  print STDERR "Group $n = $colour_of{$n}\n";
}

my @taxa = sort keys %box; # sort by --tree eventually...
my $ntaxa = scalar @taxa;
print STDERR "Number of taxa: $ntaxa\n";

my $svg = SVG->new(width=>$width+200, height=>$height);
my $row_height = int($height / @taxa);
my $row = 0;
if ($consensus) {
  label($svg, 0, 'Consensus', 'black');
  label($svg, $row_height, 'Heatmap', 'black');
  $row=2;
}

for my $t (@taxa) {

  my $Y = $row * $row_height;
  
  label($svg, $Y, $t, $colour_of{$group_of{$t}});

  for my $b (@{$box{$t}}) {
    $svg->rectangle(
      x => $b->[0] * $width / $length,
      y => $Y,
      width => $b->[1] * $width / $length,
      height => $row_height-1,
      style => { fill=>$colour_of{$b->[2]} },
    );
    
    if ($consensus) {
      $svg->rectangle(
        x => $b->[0] * $width / $length,
        y => 0,
        width => $b->[1] * $width / $length,
        height => $row_height-1,
        style => { 
          fill=>'rgb(0,0,0)',
#          'fill-opacity'=>1.0/$ntaxa,
        }
      );                                     
      $svg->rectangle(
        x => $b->[0] * $width / $length,
        y => $row_height,
        width => $b->[1] * $width / $length,
        height => $row_height-1,
        style => { 
          fill=>'rgb(0,0,0)',
          'fill-opacity'=>1.0/$ntaxa,
        }
      );                                     
    }
  }
  $row++;
}

print STDERR "Writing SVG file: $svg_fn\n";
open my $svg_fh, '>', $svg_fn;
print {$svg_fh} $svg->xmlify;
close $svg_fh;

print STDERR "Done.\n";

#----------------------------------------------------------------------

sub label {
  my($svg, $Y, $text, $colour) = @_;
  $svg->text(
    x=>$width+1, 
    y=>$Y + $row_height*0.8, 
    style=>{ 
      'font'=>'sans', 
      'font-size'=>int(0.8*$row_height),
      'font-style'=>'italic',
      'fill'=> $colour,
    },
    -cdata=>$text, 
  );
}

#----------------------------------------------------------------------

sub get_fh {
  my($fname, $desc) = @_;
  $desc ||= 'this';
  $fname or die "Please provide a $desc filename.";
  -r $fname or die "Can not open $desc file: $fname";
  print STDERR "Opening $desc: $fname\n";
  open my $fh, '<', $fname;
  return $fh;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"alignment=s",  VAR=>\$aln_fn, DEFAULT=>'', DESC=>"Alignment file"},
#    {OPT=>"tree=s",  VAR=>\$tree_fn, DEFAULT=>'', DESC=>"Tree file"},
    {OPT=>"tabular=s",  VAR=>\$tab_fn, DEFAULT=>'', DESC=>"Tabular file"},
    {OPT=>"out|svg=s",  VAR=>\$svg_fn, DEFAULT=>'', DESC=>"SVG output file"},
    {OPT=>"width=i",  VAR=>\$width, DEFAULT=>1024, DESC=>"Canvas width"},
    {OPT=>"height=i",  VAR=>\$height, DEFAULT=>600, DESC=>"Canvas height"},
#    {OPT=>"minblock=i",  VAR=>\$minblock, DEFAULT=>0, DESC=>"Boxes must be at least this wide"},
    {OPT=>"consensus!",  VAR=>\$consensus, DEFAULT=>0, DESC=>"Add consensus row"},
#    {OPT=>"border!",  VAR=>\$border, DEFAULT=>0, DESC=>"Add border to delimit genome"},
#    {OPT=>"tickmarks!",  VAR=>\$tickmarks, DEFAULT=>0, DESC=>"Add genome position tick markers"},
    {OPT=>"monochrome!",  VAR=>\$monochrome, DEFAULT=>0, DESC=>"Don't colour groups"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] -a genome.afa -t bratnextgen.tab -o pretty.svg\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__
red
green
blue
orange
brown
gray
DarkRed
LightGreen
LightBlue
gray128
Lavendar
LemonChiffon
cyan
pink

