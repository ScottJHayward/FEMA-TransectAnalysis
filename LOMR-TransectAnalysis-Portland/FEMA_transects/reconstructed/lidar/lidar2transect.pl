#!/usr/bin/env perl
use strict;
use warnings;
use lib 'c:\ourPerl';
use Mapping::UTMconvert;
use AdcircUtils::AdcGrid;
use Lidar::PointTree;
use Geometry::Interp;

#change line numbers 15, 16, and 205

#for my $nn (1..18){
# for my $nn (73,74,75,76,92,102,103,105,106,107,108,109,110,120,121,124,125,127){
for my $nn ('35.1','35.2','35.3'){

my $name = 'rc_CM-'."$nn".'XY.csv';
my $name2= 'CM-'."$nn".'XYZSTA.csv';
my $file = '../'."$name";
my $newf = "$name2";

open(my $data, '<', $file) or die "Could not open '$file' $!\n";

my @tx;
my @ty;
my @tz;
my @Z;
my @source;
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
    push @source, -999; 
	push @Z, -999; #initialize z
	print "line: $li x: $tx[$li] y: $ty[$li] z2: $Z[$li] sta: $sta[$li]\n";
  }
  $li++;  
}
$Z[3]=-999;
###################################################################

my $dir;
my $treefile;
my $sfnbin;
my $tree;
my $ds=0.00001*10;

###################################################################

$dir='F:\FEMA\LiDAR\2004-ME-Coastline_CUMBERLAND';
$treefile="$dir".'\tree.tree2';
$sfnbin="$dir".'\points.bin.fin.sfn';


#load the tree
$tree=PointTree->loadTree($treefile);
#set the superfinalized file
$tree->{SFINALIZED}=$sfnbin;


foreach my $k (0..$li-1){
   next if ($HitThisOne[$k]);
   my ($xrf,$yrf,$zrf,$idrf,$dsqrf)=$tree->getPoints_sorted($tx[$k],$ty[$k],$ds);
   print "station $sta[$k]:\n z was $Z[$k], z is $zrf->[0]\n" if ($zrf->[0]);
   next unless ($zrf->[0]);
   my $ele=$zrf->[0]*3.28084;
   $Z[$k]=$ele;
   $source[$k]='2004USGS';
   $HitThisOne[$k]=1;
   }

print "TEST 1";
##################################################################

$dir='F:\FEMA\LiDAR\2006-FEMA_CUMBERLAND';
$treefile="$dir".'\tree.tree2';
$sfnbin="$dir".'\points.bin.fin.sfn';

#load the tree
$tree=PointTree->loadTree($treefile);
#set the superfinalized file
$tree->{SFINALIZED}=$sfnbin;

foreach my $k (0..$li-1){
   next if ($HitThisOne[$k]);
   my ($xrf,$yrf,$zrf,$idrf,$dsqrf)=$tree->getPoints_sorted($tx[$k],$ty[$k],$ds);
   print "station $sta[$k]:\n z was $Z[$k], z is $zrf->[0]\n" if ($zrf->[0]);
   next unless ($zrf->[0]);
   my $ele=$zrf->[0]*3.28084;
   $Z[$k]=$ele;
   $source[$k]='2006FEMA';
   $HitThisOne[$k]=1;
}

print "TEST 2";

open(my $nf, '>', $newf);
print $nf "X,Y,Zlidar,STA,SOURCE\n";

foreach my $l (0..$li-1){
	print $nf "$tx[$l],$ty[$l],$Z[$l],$sta[$l],$source[$l]\n";
}

close $nf;

}
