cd \\SERVERME2016\Projects\2017\171.06068\Corrected_postAppeal\transect_analysis\z_correctedFiles

copyfile process_swan_output.m ../2_swan
copyfile make_WHAFIS_input.m ../3_whafis
copyfile plot_WhafisResults.m ../6_plot
copyfile READ_ME.txt ../

%re-plot the SWAN output
cd ../2_swan
delete logfiles/*_log.txt
process_swan_output

%re-run WHAFIS
cd ../3_whafis
delete whafis4/*.out 
make_WHAFIS_input

%% 
cd whafis4
% re-run the WHAFIS batch file 

%%
cd ../
shorten_output %run shorten_output file

%%
cd ../6_plot
delete *parsed.csv
%re-run the WHAFIS_out_2_kml.pl

%%
plot_WhafisResults

%% 
cd ../7_create_logs
logs2pdf

%% 
% run convert_logs.bat

%% 

move_pdf

%add a header to whafis_file:
% REVISED SEP-05-2019
% (header margin 1.5" text red, font 16)



