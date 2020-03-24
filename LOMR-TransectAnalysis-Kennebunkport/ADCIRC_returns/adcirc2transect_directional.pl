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
use Math::Trig;

# This script reads a csv file that contains lon,lat,z,station data for a wave transect 
# and then interpolates the TWL and wave conditions from the ADCIRC+SWAN model statistics
# at the transect points and then writes the transect points and interpolated
# hydrology data to another csv file to be used for further analysis.

# Input data are expected in meters (ADCIRC native), output is converted to feet. 

# the script should be configured with a list of transect id values assuming the input filename
# is formated like CM-XXX-XYZSTA.csv, where XXX is the transect id.

# this script also adjusts the TWL statistics from the ADCIRC+SWAN model to account for sea level rise
# this is necessary because YCME model used an MSL to NAVD88 conversion from VDATUM based on 
# the 1983-2001 tidal epoch to convert NACCS model boundary conditions from MSL to NAVD88
# the adjustment is a simple addition to the TWL based on the SLR rate, present year, and tidal epoch 
# configured below.

##########################################3
# start config

my @transect_list=('92','99','99-1','99-2','100','103','105','106','107','108','109','110');

my $countyName='YK-';
# ADD SEA LEVEL RISE FROM 1992 to current year

my $slrRate=1.88/1000; # long term sea level trend in meters/year from NOAA long term trend at Portland station 8418150 is 1.88 mm/yr
my $presentYear=2020;
my $epochYear=(1983+2001)/2;  # year at the midpoint of the tidal epoch
my $slrADJ=$slrRate*($presentYear-$epochYear);


#adcirc grid
my $gridFileName='fort.14';

# fort.63 file containing 100-year TWL data
my $twlFile='100-yr-TWL.63';


# end Config
#################################################################################


my $adcGrid = AdcGrid->new( $gridFileName );
my $tree = ElementQuadTree->new_from_adcGrid(-MAXELEMS=>10000,-ADCGRID=>$adcGrid);

#number of nodes
my $np=$adcGrid->getVar('NP');
my $nid=[0..$np];
my ($x,$y,$z)=$adcGrid->getNode($nid);

#read 100 year TWL
my $twl=&read63('100-yr-TWL.63');

#read 100 year PER and HS files
my @PER=();
my @HS=();
my $k=0;
foreach my $degrees (0,45,90,135,180,225,270,315){
   #Period
   my $fnameTp='100-yr-TP_'."$degrees".'.63';
   $PER[$k]=&read63($fnameTp);
   #Hs
   my $fnameHs='100-yr-HS_'."$degrees".'.63';
   $HS[$k]=&read63($fnameHs);
   #next iteration
   $k++;
}

foreach my $nn(@transect_list) {

my $name = $countyName."$nn".'XYZSTA.csv';

print "File Num: $nn \n";
print "$name \n";

my $name2 = $countyName."$nn".'XYZSTA_RETURNS.csv';

my $file = '../FEMA_transects/'."$name";
print "file: $file \n";
my $newf = $name2; #change directory if needed

open(my $data, '<', $file) or die "Could not open '$file' $!\n";

my @tx=();
my @ty=();
my @tz=();
my @Z=();
my @sta=();
my @HitThisOne=();  # array to tell us if we already have data for this node
my $li=-1;

while (my $line = <$data>) {
  chomp $line;
  if ($li > -1){
    my @fields = split "," , $line;
    $tx[$li]=$fields[0];
    $ty[$li]=$fields[1];
    $tz[$li]=$fields[2];
	$sta[$li]=$fields[3];
	push @Z, -999; #initialize z
  }
  $li++;  
}
$Z[3]=-999;

# find the heading of the transect
my $v_x = $tx[$#tx]-$tx[1];
my $v_y = $ty[$#ty]-$ty[1];
my $pi = 3.14159;
my $heading = 180*atan2($v_y,$v_x)/$pi;
if ($heading < 0) {
    $heading=$heading+360;
}



#decide which directional bin
my $per=$PER[0]; 
my $hs=$HS[0];
my $dir=0;

if ($heading > 22.5){
$per=$PER[1];
$hs=$HS[1];
$dir=45;
}
 if ($heading > 67.5){
$per=$PER[2];
$hs=$HS[2];
$dir=90;
}
if ($heading > 112.5){
$per=$PER[3];
$hs=$HS[3];
$dir=135;
}
if ($heading > 157.5){
$per=$PER[4];
$hs=$HS[4];
$dir=180;
}
if ($heading > 202.5){
$per=$PER[5];
$hs=$HS[5];
$dir=225;
} 
if ($heading > 247.5){
$per=$PER[6];
$hs=$HS[6];
$dir=270;
}
if ($heading > 292.5){
$per=$PER[7];
$hs=$HS[7];
$dir=315;
} 
if ($heading > 337.5){
$per=$PER[0];
$hs=$HS[0];
$dir=0;
} 


open(my $nf, '>', $newf);
print $nf "X,Y,Z,STA,Zadcirc,HSadcirc,TWLadcirc,PERadcirc,TransHeading,PerHeading\n";

foreach my $l (0.. $li-1){
    my $lon=$tx[$l];
    my $lat=$ty[$l];
    my $interpolant = $tree->getInterpolant(-XX=>$lon,-YY=>$lat);
    #my $val = $tree->getZvalue(-ZDATA=>$z,-XX=>$lon,-YY=>$lat);
    my $val =-99999;
    $val = $tree->interpValue(-ZDATA=>$z,-INTERPOLANT=>$interpolant) if defined $interpolant;
    my $zfeet = -3.28083*$val;
    #my $val2 = $tree->getZvalue(-ZDATA=>\@hs,-XX=>$lon,-YY=>$lat);
    my $val2=-99999;
    $val2 = $tree->interpValue(-ZDATA=>$hs,-INTERPOLANT=>$interpolant) if defined $interpolant;
    my $HSfeet = 3.28083*$val2;
    #my $val3 = $tree->getZvalue(-ZDATA=>\@twl,-XX=>$lon,-YY=>$lat);
    my $val3=-99999;
    $val3 = $tree->interpValue(-ZDATA=>$twl,-INTERPOLANT=>$interpolant) if defined $interpolant;
    my $TWLfeet = 3.28083*($val3+$slrADJ);
    #my $period = $tree->getZvalue(-ZDATA=>\@per,-XX=>$lon,-YY=>$lat);
    my $period=-99999;
    $period = $tree->interpValue(-ZDATA=>$per,-INTERPOLANT=>$interpolant) if defined $interpolant;
    print "Transect: $nn , Heading: $heading \n  line: $l x: $tx[$l] y: $ty[$l] z: $tz[$l] sta: $sta[$l] \n  z(adc): $zfeet hs: $HSfeet twl: $TWLfeet per: $period  \n";
    print $nf "$tx[$l],$ty[$l],$tz[$l],$sta[$l],$zfeet,$HSfeet,$TWLfeet,$period,$heading,$dir\n";
}

close $nf;

}



sub read63{

   my $filename = shift; 
   open(my $fh, $filename)
    or die "Could not open file '$filename' $!";
   # skip first three lines
   <$fh>;<$fh>;<$fh>;
   my @TWL;
   while (my $row = <$fh>) {
     chomp $row;
     $row =~ s/^\s+//; # remove leading whitespace 
     my @data=split " " , $row;
     $TWL[$data[0]]=$data[1];
   }

   return \@TWL;
}

