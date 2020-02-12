#!/usr/bin/env perl
use strict;
use warnings;
# use Text::CSV;
use lib 'c:\ourPerl';
use Mapping::UTMconvert;
use AdcircUtils::AdcGrid;
use AdcircUtils::ElementQuadTree;
use Lidar::PointTree;
use Geometry::Interp;

#change line numbers 15, 16, and 205

#adcirc grid
my $gridFileName='fort.14';
my $adcGrid = AdcGrid->new( $gridFileName );
my $tree = ElementQuadTree->new_from_adcGrid(-MAXELEMS=>10000,-ADCGRID=>$adcGrid);

#number of nodes
my $np=$adcGrid->getVar('NP');
my $nid=[0..$np];
my ($x,$y,$z)=$adcGrid->getNode($nid);


#read 100 year Hs
my $filename = 'HS_100.txt'; 
open(my $fh, $filename)
  or die "Could not open file '$filename' $!";
my @hs=0;
my $line=0;
while (my $row = <$fh>) {
  chomp $row;
  my @data=split " " , $row;
  $hs[$line]=$data[1];
  print "$hs[$line]\n";
  $line++;
}

#read 100 year TWL
my $filename2 = 'TWL_100.txt'; 
open(my $fh2, $filename2)
  or die "Could not open file '$filename' $!";
my @twl=0;
my $line=0;
while (my $row = <$fh2>) {
  chomp $row;
  my @data=split " " , $row;
  $twl[$line]=$data[1];
  print "$twl[$line]\n";
  $line++;
}


# for my $nn (73,74,75,76,92,102,103,105,106,107,108,109,110,120,121,124,125,127){
for my $nn ('35.3'){

my $file = '../rc_CM-'."$nn".'XY.csv';
my $newf = 'CM-'."$nn".'XYZSTA.csv';

open(my $data, '<', $file) or die "Could not open '$file' $!\n";

my @tx;
my @ty;
my @tz;
my @Z;
my @sta;
my @HitThisOne;  # array to tell us if we already have data for this node
my $li=-1;

while (my $line = <$data>) {
  chomp $line;
  if ($li > -1){
    my @fields = split "," , $line;
    $tx[$li]=$fields[0];
    $ty[$li]=$fields[1];
	$sta[$li]=$fields[2];
	push @Z, -999; #initialize z
  }
  $li++;  
}
$Z[3]=-999;


open(my $nf, '>', $newf);
print $nf "X,Y,Zadcirc,STA,source\n";

foreach my $l (0.. $li-1){
    my $lon=$tx[$l];
    my $lat=$ty[$l];
    my $val = $tree->getZvalue(-ZDATA=>$z,-XX=>$lon,-YY=>$lat);
    my $zfeet = -3.28083*$val;
    print "line: $l x: $tx[$l] y: $ty[$l] z: $zfeet \n";
    print $nf "$tx[$l],$ty[$l],$zfeet,$sta[$l],interpolated ADCIRC\n";
}

close $nf;

}
