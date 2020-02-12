2017 FEMA Appeal: Transect-Based Analyses
___________________________________________________________________________


This document describes how to set up and run the transect based hydrodynamic 
analyses for the York and Cumberland County FEMA appeals, prepared by Scott 
Hayward, Nathan Dill, and Tom Nielson of Ransom Conlsulting inc. 

nld-20190912


___________________________________________________________________________
Transect Preparation

The data in the preliminary flood insurance study provided by FEMA were not 
consistent on a town-by-town basis. The most notable inconsistency being the 
resolution and formatting of the provided 1-dimensional elevation profiles. 
Before configuring or proceeding with any transect-based calculations the engineer
should asses the provided profile, and determine whether the horizontal spatial 
reference is sufficient to be used for hydrodynamic analyses. The Directory titled
'FEMA_transects' should contain all transects used for this appeal. 

Here are the most common scenarios at this stage: 

1. FEMA has provided the full elevation profile, along with horizontal coordinates.
The transect appears to be representative of the shoreline. 

2. The full transect was not provided, but a lower resolution transect (from CHAMP) 
was provided, containing horizontal reference data. 

3. A full resolution profile was provided, but horizontal reference was not available. 

4. Sufficient data could not be found in the preliminary flood insurance study,
or anew transect was created. 

For scenario 1, a CSV file for each of FEMA's transsects should be placed in 
the 'FEMA_transects' directory. The columns should be:
X(degrees longitude), Y(degrees latitude), Z(ft), Station(ft)

For the remaining options, the user should use either the matlab scripts in 
the  'WHAFIS_only' or the 'nospatial' directories to create a new profile where
elevation data should be interpolated. Once the profiles have been made, the 
transects will be saved in the 'reconstructed' directory. Perl scripts in the 
LIDAR and BATHY directories are used to interpolate the most recently avaiable 
elevation data onto the transect. A matlab script called ADC_LiDAR is used 
to combine the datasets, and write the profiles into the correct format in 
the FEMA_transects directory. 

Finally, configure the perl script in ADCIRC_returns for the transects used 
in the appeal. A new .csv file will be created for each transect, and will 
contain the statistics from the hydrologic analysis interpolated onto each 
point along the elevation profile. 

___________________________________________________________________________
1: User input
Configure the matlab script in this directory to run on desired transects. 
There is also an array defined below the transect names that allows this script
to be run on specific transects without changing the toe/top locations. 

This script will plot the adcirc returns along the elevation profile, and prompt the user
to select the starting point for WHAFIS/SWAN and then the toe/top locations for eah transect. 
It is important that the user selects a WHAFIS/SWAN start that contains wave data. If there are 
no wave data returns along the transect, it can be added later by opening up transectdata.csv 
and manually changing the results in the appropriate columns befire running SWAN and WHAFIS. 

___________________________________________________________________________
2: SWAN

- In Matlab change directory to the 2_swan directory

- If necessary, delete unwanted files in the logfiles, swanfiles, and gridfiles directory

- Run the make_swan_input.m script.

- In Windows Explorer change to the swanfiles directory and run runswan.bat

- In Matlab run the process_swan_outpt.m script

- Check the plots and logs in the logfiles directory

___________________________________________________________________________
3: WHAFIS

- Change to the 3_whafis directory

- If necessary, delete files in the logfiles and whafis4 directory, 
  but do not delete the MG.dat and WHAFIS4 executable

- In Windows, change to the whafis4 directory and run runWHAFIS.bat

- In Matlab, run shorten_output.m (this removes blank lines from WHAFIS output files) 

- Optional: go to step #6 to plot WHAFIS results 
___________________________________________________________________________
4: TAW

- Change to the 4_taw directory

- In Matlab run the TAW_iterative_writer.m script

- Open the Taw_iterative.m script in a text editor and copy and paste into the
  Matlab command window. If you do not copy past the diary wont be written correctly

- In Matlab, run the shorten_Diary.m script to remove blank lines from the diary file

- In Matlab, run the TAW_logfiles.m script to create the log files

- Review TAW diary, logs, and plots,  if deisired adjust toe/top station in ../data/transectdata.xls
  and go to step #1

____________________________________________________________________________
5: RUNUP2

- Change to the 5_runup2 directory
 
- If necessary, delete files in logfiles and RUNUP2_files directories,
  but do not delete RUNUP2.exe executable

- In Matlab, run make_runup2_input.m

- Change to the RUNUP2_files directory and run run.bat 
  (note, RUNUP2.exe will generally not run on modern 64-bit machines/OS, 
  we've found that the MS DOS emulator DOSBOX is a conveient way to run
  this ancient program www.dosbox.com)

- In Matlab, run the process_runup2_output.m script.  This also runs
  the ACES_Beach_Runup.m script and reports the output in the logfile

___________________________________________________________________________
6: Plot

- Change to the 6_plot directory

- Run the WHAFIS_out_2_KML.pl script
  this script requres that you have the ourPerl library available at c:\ourPerl
  ourPerl can be obtained from github http://github.com/natedill/ourPerl.git
  it is also provided on the enclosed DVD 

- In Matlab, run the plot_WhafisResults.m script

- In Matlab, run the runupZones.m script

- The above scripts generate kml files that allow you to view the WHAFIS
  and RUNUP2/TAW results in Google Earth

____________________________________________________________________________
7: Calculation Logs

- Change to the 7_create_logs directory

- In Matlab, run the logs2pdf.m script

- In Windows run the convert_logs.bat batch file, this requires the txt2pdf.win95.exe program

- In Matlab, run the move_pdf.m script. This will copy all the calculation records into a 
  folder for each transect evaluated.  Input and results from all calculations will also 
  be sumarized in the data/transectdata.xls MS Excel file











