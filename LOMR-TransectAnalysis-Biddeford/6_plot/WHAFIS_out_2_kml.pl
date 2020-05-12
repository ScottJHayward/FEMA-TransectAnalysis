#!/usr/bin/env perl
use strict;
use warnings;
use lib 'c:\ourPerl';
use Mapping::UTMconvert;

my $whafisdir='../3_whafis/whafis4/';

opendir( DIR, $whafisdir);

my @FILES=readdir(DIR);

foreach my $Tname (@FILES){

 next unless $Tname =~ m/\.out$/;
 my $whafisname=$Tname;
 $whafisname =~ s/\.out$//;
 my $whafisOutFile="$whafisdir"."$whafisname".'.out';
 my $kmlFile="$whafisname".'.kml';
 my $csvFile="$whafisname"."_parsed.csv";
 my $spZone='19 T';

 # get the start and end points (assume a PS card at the end) e.g.
 #                         PS#  1   start(2645889.240,465500.797) end(2646212.790,469238.578)       
 # read the WHAFIS data
 open FILE, "<$whafisOutFile"; 

 my @start;
 my @end;
 while (<FILE>){
   chomp;
   if ($_ =~ m/end\(([0-9,\.]+)\)/i) {
     @end=split(/,/,$1);
   }
   if ($_ =~ m/start\(([0-9,\.]+)\)/i) {
     @start=split(/,/,$1);
   }
 }
 close(FILE);

 print " start = @start; end = @end\n";

   my @start_l=UTMconvert::utm2deg($start[0],$start[1],$spZone);
   my @end_l=UTMconvert::utm2deg($end[0],$end[1],$spZone);

 my $dx=$end[0]-$start[0]; 
 my $dy=$end[1]-$start[1]; 
 my $len=($dx**2 + $dy**2)**0.5; 
 my $dldx=$len/$dx;
 my $dldy=$len/$dy;

 # get the part 1 data and store in CARDS array of hash references
 open FILE, "<$whafisOutFile";
 my @CARDS;
 while (<FILE>){
   my %params=();
   chomp;
   $_ =~ s/^\s+//;
   my $type=substr($_,0,2);
   
   $params{type}=$type;
   
   my $sta;
   my $elev;   
   my $owFetch;
   my $SWL_10;
   my $SWL_100;
   my $Hs;
   my $Tp;

   if     (uc($type) eq 'IE'){  # IE CARD
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8);
       #$owFetch=substr ($_,16,8);
       #$SWL_10=substr ($_,24,8);
       $SWL_100=$data[4];#substr ($_,32,8);
       #$Hs=substr ($_,40,8);
       #$Tp=substr ($_,48,8);
   }elsif (uc($type) eq 'AS' ){   # AS CARD
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8);       
       #$sta=substr ($_,2,6);
       #$elev=substr ($_,8,8);
       #$SWL_10=substr ($_,16,8);
       $SWL_100=$data[3];#substr ($_,24,8);
   }elsif (uc($type) eq 'BU' ){  # BU CARD
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8); 
      # $sta=substr ($_,2,6);
      # $elev=substr ($_,8,8);
      # $SWL_10=substr ($_,32,8);
       $SWL_100=$data[5];#substr ($_,40,8);
   }elsif (uc($type) eq 'DU' ){ 
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8); 
       #$sta=substr ($_,2,6);
      # $elev=substr ($_,8,8);
      # $SWL_10=substr ($_,24,8);
       $SWL_100=$data[4];#substr ($_,32,8);
   }elsif (uc($type) eq 'IF' ){ 
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8); 
      # $sta=substr ($_,2,6);
      # $elev=substr ($_,8,8);
      # $SWL_10=substr ($_,16,8);
       $SWL_100=$data[3];#substr ($_,24,8);
   }elsif (uc($type) eq 'OF' ){ 
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8); 
      # $sta=substr ($_,2,6);
      # $elev=substr ($_,8,8);
      # $SWL_10=substr ($_,16,8);
       $SWL_100=$data[3];#substr ($_,24,8);
   }elsif (uc($type) eq 'VE' ){ 
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8); 
       #$sta=substr ($_,2,6);
       #$elev=substr ($_,8,8);
       #$SWL_10=substr ($_,48,8);
       $SWL_100=$data[7];#substr ($_,56,8);
   }elsif (uc($type) eq 'VH' ){ 
       my @data=split(/\s+/,$_);
       shift @data;
       $sta=$data[0];#substr ($_,2,6);
       $elev=$data[1];#substr ($_,8,8); 
       #$sta=substr ($_,2,6);
       #$elev=substr ($_,8,8);
       #$SWL_10=substr ($_,48,8);
       $SWL_100=$data[7];#substr ($_,56,8);
   }elsif (uc($type) eq 'MG ' ){ 
       next;
   }elsif (uc($type) eq 'ET' ){ 
       last;
   }else{
       next;
   
   }
   
   #convert to m
   $sta=$sta/3.280833333;
   
   $params{station}=$sta;
   $params{elevation}=$elev;
   $params{SWL_100}=$SWL_100;
  # $params{SWL_10}=$SWL_10;
   
   my $x=($start[0] + $sta/$dldx);
   my $y=($start[1] + $sta/$dldy);


   # get lon and lat
   my ($lon,$lat)=UTMconvert::utm2deg($x,$y,$spZone);
   # $lon=-$lon;
   $params{lon}=$lon;
   $params{lat}=$lat;
   $params{x}=$x;
   $params{y}=$y;


   push @CARDS, \%params;

   print "$sta, $x, $y, $elev, $SWL_100\n";

 }
 close(FILE);

   

 # parse part 2
 my @P2;

 # open FILE, "<$whafisOutFile";
 # my @CARDS2;
 # while (<FILE>){
   # chomp;
   # if ($_ =~ m/PART2/){
      # <FILE>; # skip five lines
      # <FILE>;
      # <FILE>;
      # <FILE>;
      # <FILE>;
      # foreach my $c (@CARDS){
         # my %params=%{$c};
         # my $line=<FILE>;
         # chomp ($line);
         # push @P2, $line;
         # $line =~ s/^\s+//;
         # my ($cd,$sta,$Hc,$Tp,$Crest)=split(/\s+/,$line);
         # # check sta to see if it matches
         # $params{Hc}=$Hc;
         # $params{Tp}=$Tp;
         # $params{Crest}=$Crest;
         # <FILE>;  # because of the useless blank lines in WHAFIS part 2 output
         # push @CARDS2, \%params;

      # }
   # }
 # }



 open FILE, "<$whafisOutFile";
 my @CARDS2;
 while (<FILE>){
   chomp;
   if ($_ =~ m/PART2/){
      <FILE>; # skip five lines
      <FILE>;
      <FILE>;
      <FILE>;
      <FILE>;
      foreach my $c (@CARDS){
         my %params=%{$c};
         my $line;
         print "$params{type}\n";
         while (<FILE>){
            $line=$_;
            print "$line";
            if ($line =~ m/$params{type}/){
               last;
            } else {
               next;
            }
         }
         chomp ($line);
         push @P2, $line;
         $line =~ s/^\s+//;
         my ($cd,$sta,$Hc,$Tp,$Crest)=split(/\s+/,$line);
         # check sta to see if it matches
         $params{Hc}=$Hc;
         $params{Tp}=$Tp;
         $params{Crest}=$Crest;
         <FILE>;  # because of the useless blank lines in WHAFIS part 2 output
         push @CARDS2, \%params;
      }
   }
 }
 close (FILE);
     

 #my $tt=$CARDS2[10];
 #my %hh=%{$tt};

 #print "jjjjjjjj:::::$hh{type} $hh{station} $hh{SWL_100} $hh{Hc}\n";


 open KML, ">$kmlFile"; 

 print KML"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n<Document>\n";
 print KML "   <name>$kmlFile</name>\n";
 print KML "   <Style id=\"Tstyle\">\n";
 print KML "     <IconStyle>\n";
 print KML "        <scale>0.6</scale>\n";
 # print KML "            <Icon>\n";
 # print KML "               <href>dot.png</href>\n";
 # print KML "            </Icon>\n";
 print KML "     </IconStyle>\n";
 print KML "     <LineStyle><width>3</width></LineStyle>\n";
 print KML "    </Style>\n";
 print KML " <Folder>\n";
 print KML "  <name>$Tname</name>\n";
 print KML "  <Placemark>\n";
 print KML "    <name>$Tname</name>\n";
 print KML "    <description>\n";
 print KML "    PART 2 WHAFIS output: \n\n";
 foreach my $l (@P2){
     print KML "$l\n";
 }
 print KML "    </description>\n";
 print KML '    <styleUrl>#Tstyle</styleUrl>';
 print KML "\n";
 print KML "    <LineString>\n";
 print KML "    <tessellate>1</tessellate>\n";
 print KML "    <coordinates>\n";
 print KML "      $start_l[0],$start_l[1],0 $end_l[0],$end_l[1],0\n";	#			-93.3089692412841,29.7540426813781,0 -93.30401219050678,29.83080666869899,0 
 print KML "    </coordinates>\n";
 print KML "    </LineString>\n";
 print KML "  </Placemark>\n";
 print KML "</Folder>\n";

 print KML " <Folder><name>WHAFIS CARDS</name>\n";
 foreach my $c (@CARDS2){
     

     my  %cc=%{$c};


     my $lon=$cc{lon};
     my $lat=$cc{lat};
     my $hc=$cc{Hc};
     my $crest=$cc{Crest};
     my $bfe=int($crest+0.5);
     my $zone='AE';
     $zone='VE' if ($hc >=3);
     next if $cc{type} eq 'AS';
     $cc{zone}=$zone;
     $cc{bfe}=$bfe;




     print KML "  <Placemark>\n";
     #print KML "    <name>$cc{type} station $cc{station} ft</name>\n";
     print KML "    <name>$zone $bfe</name>\n";
     print KML "    <description>\n";
     foreach my $key (keys %cc){
         print KML "$key :: $cc{$key};";
     }
     print KML "    </description>\n";
     print KML '    <styleUrl>#Tstyle</styleUrl>';
     # print KML "    <visibility>0</visibility>\n";
     print KML "\n";
     print KML "    <Point>\n";
     print KML "    <coordinates>\n";
     print KML "      $lon,$lat,$cc{elevation}\n";	#			-93.3089692412841,29.7540426813781,0 -93.30401219050678,29.83080666869899,0 
     print KML "    </coordinates>\n";
     print KML "    </Point>\n";
     print KML "  </Placemark>\n";
 }
 print KML "</Folder>\n";


 print KML "</Document>\n";
 print KML "</kml>\n";
 close(KML);



 open CSV, ">$csvFile";
 print CSV "type,station-ft,x-ft,y-ft,lon,lat,elev-ft,SWL-ft,Hc-ft,CrestElev-ft,Tp\n";

 foreach my $c(@CARDS2){
     my  %cc=%{$c};
     my $sta=$cc{station}*3.280833333;

     print CSV "$cc{type},$sta,$cc{x},$cc{y},$cc{lon},$cc{lat},$cc{elevation},$cc{SWL_100},$cc{Hc},$cc{Crest},$cc{Tp}\n";
 }
 close (CSV);





       


} # end loop over files
