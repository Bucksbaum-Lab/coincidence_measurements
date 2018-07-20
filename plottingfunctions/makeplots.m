clear all
close all

f = momentumPlotsFunctions();

%put in data info
folder = 'G:\2018_05_11\combined\';
filename = 'acetylene_1ps_6p4en10_torr_800_1300_266';
delay = '2000';
intensity = 'low';

shutterMarkers = [string('_closed_'), string('_open_')];
shutterTitle = [string('control off'), string('control on')];

delayMarkers = [string('_p'), string('_n')];
delayTitle = [string(['+' delay ' fs']), string(['-' delay ' fs'])];

%move to center of mass frame
momsumlimit = 2;

[plotting_data] = ...
    getmomentum(folder, filename, intensity, delay, momsumlimit, shutterMarkers, delayMarkers);

%update the energy
plotting_data = resetEV(plotting_data);

%project the hydrogen momentum onto the CH+ and C+ ions
plotting_data = anglesAndProjections(plotting_data);

%make polar plots
numMom = 10;
maxMom = 8;
numAngle = 5;
%%
plotting_data = polarPlots(plotting_data, numMom, maxMom, numAngle, shutterTitle, delayTitle, intensity);
plotting_data = polarPlotsNormalized(plotting_data, numMom, maxMom, numAngle, shutterTitle, delayTitle, intensity);
%%
compareShutterAngle(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, intensity)
compareShutterLowAnulus(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, numMom, maxMom, intensity)
compareShutterHighAnulus(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, numMom, maxMom, intensity)
%%
compareDelayAngle(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, intensity)
compareDelayLowAnulus(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, numMom, maxMom, intensity)
compareDelayHighAnulus(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, numMom, maxMom, intensity)
%%
%pcolor(plotting_data(1).Xpolar, plotting_data(1).Ypolar, plotting_data(1).polarDist)
%%
for ii = 1:4
    [plotting_data(ii).clusters, plotting_data(ii).centroid] = kmeans([plotting_data(ii).parallel_proj,plotting_data(ii).perpendicular_proj],4);
    figure();gscatter(plotting_data(ii).parallel_proj,plotting_data(ii).perpendicular_proj,plotting_data(ii).clusters)
    title([plotting_data(ii).shutter, plotting_data(ii).delay])
end