function varargout = main(varargin)

%{
|        |                  |         |
|        |                  |         |
|     <--|   Flight Tube    |         | MCP
|      s |                  |         |

V1 <---> G <--------------> G <-----> VM
    A    r        L         r    C
         o                  o
         u                  u
         n                  n
         d                  d
%}

% MAIN MATLAB code for main.fig
%    MAIN, by itself, creates a new MAIN or raises the existing
%    singleton*.
%
%    H = MAIN returns the handle to a new MAIN or the handle to
%    the existing singleton*.
%
%    MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%    function named CALLBACK in MAIN.M with the given input arguments.
%
%    MAIN('Property','Value',...) creates a new MAIN or raises the
%    existing singleton*.  Starting from the left, property value pairs are
%    applied to the GUI before main_OpeningFcn gets called.  An
%    unrecognized property name or invalid value makes property application
%    stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%    *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%    instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 11-Jul-2018 17:23:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    try
        gui_mainfcn(gui_State, varargin{:});
    catch ME
        errordlg(ME.message,'Error');
        ME.rethrow
    end
end
% End initialization code - DO NOT EDIT



% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

set(handles.chunksStartNum, 'string', '0');

% Update handles structure
guidata(hObject, handles);


function file = loadIfExists(filename, length)
    
    if(~exist('length', 'var')) length = 2; end
    
    if(exist(filename, 'file') == 2)
        file = load(filename);
    else
        warning([filename ' does not exist. Loaded "0".'])
        warndlg([filename ' does not exist. Loaded "0".'])
        file = zeros(length, 1);
    end

    
% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

loadSave = get(handles.loadSave, 'value');

%set the error bar to no errors
set(handles.err, 'string', 'please select a file, raw or pre-loaded');

%check that the file exist
[file, handles.path] = uigetfile({'*.txt;*.mat'});
set(handles.err, 'string', 'no errors');
set(handles.loadfile, 'string', 'loading');
drawnow

filename = [handles.path file];

[~, handles.datafile, ext] = fileparts(filename);

if strcmp(ext, '.txt')
    [handles, delayInfo, polarInfo] = loadTxtData(handles, filename);
elseif strcmp(ext, '.mat')
    [handles, delayInfo, polarInfo] = loadMatData(handles, filename);
else
    set(handles.err, 'string', 'you did not select a .mat or .txt file');
    error('you did not select a .mat or .txt file')
end
    
if loadSave
    saveLoadedData(handles, delayInfo, polarInfo)
end

%save any values saved to handles
guidata(hObject, handles);

clearvars -except handles

set(handles.loadfile, 'string', 'load');

function [handles, delayInfo, polarInfo] = loadTxtData(handles, filename)
                                                  
breakUpData = get(handles.breakUpData, 'value');
numChunks = get(handles.numEventTags, 'string');
numChunks = str2double(numChunks);

%estimate how long it will take to load file
estimatedTime = 0;

for ii = 1:10
    tic
    testLoad = load([cd '\10000points.dat']);
    estimatedTime = estimatedTime+toc;
end

fileInfo = dir(filename);
fileInfoTest = dir([cd '\10000points.dat']);
estimatedTime = estimatedTime/fileInfoTest.bytes/10*fileInfo.bytes;
clear testLoad;

set(handles.err, 'string', ['estimated completed load time: '...
    datestr(datetime('now')+seconds(estimatedTime))]);
drawnow

%load the file and the other info required
if exist([handles.path '\DelayInfo.txt'], 'file') == 2
    delayInfo = load([handles.path '\DelayInfo.txt']);
    x = cell(size((delayInfo(1,:))',1)+1,1);
    for nn = 1:size((delayInfo(1,:))',1)
        x(nn+1) = {num2str(delayInfo(1,nn))};
    end
    x(1) = {'all'};

    set(handles.delayChoice, 'string', x);
else
    set(handles.delayChoice, 'string', {'all'})
    delayInfo = 0;
end

handles.delayInfo = delayInfo;

if exist([handles.path '\polarizationInfo.txt'], 'file') == 2
    polarInfo = load([handles.path '\polarizationInfo.txt']);
    x = cell(size((polarInfo(1,:))',1)+1,1);
    for nn = 1:size((polarInfo(1,:))',1)
        x(nn+1) = {num2str(polarInfo(1,nn))};
    end
    x(1) = {'all'};

    set(handles.polarChoice, 'string', x);
else
    set(handles.polarChoice, 'string', {'all'})
    polarInfo = 0;
end

handles.polarInfo = polarInfo;

if exist([handles.path '\parameterInfo.txt'], 'file') == 2
    paramInfo = readtable([handles.path '\parameterInfo.txt']);
    x = cell(size(paramInfo(2,:),2)+1,1);
    x(2:end) = table2cell(paramInfo(2,:));
    x(1) = {'all'};
    x = lower(x);

    set(handles.paramChoice, 'string', x);
else
    set(handles.paramChoice, 'string', {'all'})
    paramInfo = 0;
end

handles.paramInfo = paramInfo;

eventtags = loadIfExists([handles.path '\eventtags.txt']);
shutterstatus = loadIfExists([handles.path '\shutterclosed.txt'], size(eventtags, 1));
delaystatus = loadIfExists([handles.path '\Delay.txt'], size(eventtags, 1));
polarizationstatus = loadIfExists([handles.path '\Polarization.txt'], size(eventtags, 1));
paramstatus = loadIfExists([handles.path '\ParamIndex.txt'], size(eventtags, 1));
intensitystatus = loadIfExists([handles.path '\Intensity.txt'], size(eventtags,1));

if breakUpData

    chunksStartNum = str2double(get(handles.chunksStartNum, 'string'));

    timee = clock;
    loadInfo = [date, '-', num2str(timee(4)), '-', num2str(timee(5)), '-', num2str(floor(timee(6))), '_'];

    rawDataStore = tabularTextDatastore(filename);
    EventTagn = [eventtags; 10^10];

    loadSave = false;

    for mm = 1:floor((length(EventTagn)-1)/numChunks)

        if hasdata(rawDataStore)

            rawDataStore.ReadSize = EventTagn(numChunks*mm+1)-EventTagn(numChunks*mm-(numChunks-1));
            rawdata = table2array(read(rawDataStore));

            %initialize the arrays for holding ion info
            [numshots,maxions] = size(rawdata);
            maxions = (maxions-2)/3;

            desiredSize = rawDataStore.ReadSize;

            while (numshots ~= desiredSize) && hasdata(rawDataStore)

                rawDataStore.ReadSize = rawDataStore.ReadSize - numshots;
                rawdata = [rawdata; table2array(read(rawDataStore))];

                [numshots, ~] = size(rawdata);

            end

            eventtags = [EventTagn(((mm-1)*numChunks+1):mm*numChunks)-EventTagn((mm-1)*numChunks+1); numshots];
            
            ions_x = zeros(numshots, maxions);
            ions_y = zeros(numshots, maxions);
            ions_tof = zeros(numshots, maxions);

            %save the info into respective arrays
            handles.numHitsRaw = rawdata(:, 2);

            for nn = 1:maxions
                ions_x(:, nn) = rawdata(:, 3*nn);
                ions_y(:, nn) = rawdata(:, 3*nn+1);
                ions_tof(:, nn) = rawdata(:, 3*nn+2);
            end

            handles.ions_x = ions_x;
            handles.ions_y = ions_y;
            handles.ions_tof = ions_tof;

            shutterstatusRaw = zeros(numshots, 1);
            delaystatusRaw = zeros(numshots, 1);
            polarizationstatusRaw = zeros(numshots, 1);
            paramstatusRaw = zeros(numshots, 1);
            intensitystatusRaw = zeros(numshots, 1);

            for nn = 1:length(eventtags)-1

                shutterstatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
                    repmat(shutterstatus(((mm-1)*numChunks+1)+nn-1), (eventtags(nn+1)-eventtags(nn)), 1);

                delaystatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
                    repmat(delaystatus(((mm-1)*numChunks+1)+nn-1), (eventtags(nn+1)-eventtags(nn)), 1);

                polarizationstatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
                    repmat(polarizationstatus(((mm-1)*numChunks+1)+nn-1), (eventtags(nn+1)-eventtags(nn)), 1);

                paramstatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
                    repmat(paramstatus(((mm-1)*numChunks+1)+nn-1), (eventtags(nn+1)-eventtags(nn)), 1);

                intensitystatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
                    repmat(intensitystatus(((mm-1)*numChunks+1)+nn-1), (eventtags(nn+1)-eventtags(nn)), 1);

            end

            handles.shutterstatusRaw = shutterstatusRaw;
            handles.delaystatusRaw = delaystatusRaw;
            handles.polarizationstatusRaw = polarizationstatusRaw;
            handles.paramstatusRaw = paramstatusRaw;
            handles.intensitystatusRaw = intensitystatusRaw;

            loaded_data = struct('datafile', handles.datafile, 'path', handles.path,...
                'numHitsRaw', handles.numHitsRaw, 'ions_x', handles.ions_x,...
                'ions_y', handles.ions_y, 'ions_tof', handles.ions_tof,...
                'shutterstatusRaw', handles.shutterstatusRaw,...
                'delaystatusRaw', handles.delaystatusRaw,...
                'polarizationstatusRaw', handles.polarizationstatusRaw,...
                'paramstatusRaw', handles.paramstatusRaw,...
                'intensitystatusRaw', handles.intensitystatusRaw,...
                'shotsStartZero', EventTagn((mm-1)*numChunks+1),...
                'delayInfo', delayInfo, 'polarInfo', polarInfo, 'paramInfo', paramInfo);

            if ~exist([handles.path '\analysis\loadedData'], 'dir')
                mkdir([handles.path '\analysis\loadedData']);
                pause(1)
            end

            filename = [handles.datafile, '-loaded_data-', loadInfo, num2str(chunksStartNum+mm), '.mat'];
            filename = fullfile([handles.path '\analysis\loadedData'], filename);

            save(filename, 'loaded_data', '-v7.3');

        end
    end

else

    rawdata = sparse(load(filename));

    %initialize the arrays for holding ion info
    [numshots,maxions] = size(rawdata);
    eventtags = [eventtags; numshots];
    maxions = (maxions-2)/3;

    ions_x = zeros(numshots, maxions);
    ions_y = zeros(numshots, maxions);
    ions_tof = sparse(numshots, maxions);

    %save the info into respective arrays
    handles.numHitsRaw = rawdata(:, 2);

    for nn = 1:maxions
        ions_x(:, nn) = rawdata(:, 3*nn);
        ions_y(:, nn) = rawdata(:, 3*nn+1);
        ions_tof(:, nn) = rawdata(:, 3*nn+2);
    end

    handles.ions_x = ions_x;
    handles.ions_y = ions_y;
    handles.ions_tof = ions_tof;

    shutterstatusRaw = nan(eventtags(length(eventtags)), 1);
    delaystatusRaw = nan(eventtags(length(eventtags)), 1);
    polarizationstatusRaw = nan(eventtags(length(eventtags)), 1);
    paramstatusRaw = nan(eventtags(length(eventtags)), 1);
    intensitystatusRaw = nan(eventtags(length(eventtags)), 1);

    for nn = 1:length(eventtags)-1
        shutterstatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
            repmat(shutterstatus(nn), (eventtags(nn+1)-eventtags(nn)), 1);
        delaystatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
            repmat(delaystatus(nn), (eventtags(nn+1)-eventtags(nn)), 1);
        polarizationstatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
            repmat(polarizationstatus(nn), (eventtags(nn+1)-eventtags(nn)), 1);
        paramstatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
            repmat(polarizationstatus(nn), (eventtags(nn+1)-eventtags(nn)), 1);
        intensitystatusRaw((eventtags(nn) + 1):eventtags(nn + 1)) = ...
            repmat(intensitystatus(nn), (eventtags(nn+1)-eventtags(nn)), 1);

    end

    handles.shutterstatusRaw = sparse(shutterstatusRaw);
    handles.delaystatusRaw = sparse(delaystatusRaw);
    handles.polarizationstatusRaw = sparse(polarizationstatusRaw);
    handles.paramstatusRaw = sparse(paramstatusRaw);
    handles.intensitystatusRaw = intensitystatusRaw;

end

function [handles, delayInfo, polarInfo] = loadMatData(handles, filename)
    load(filename)
    
    handles.datafile = loaded_data.datafile;
    handles.path = loaded_data.path;
    handles.numHitsRaw = loaded_data.numHitsRaw;
    handles.ions_x = loaded_data.ions_x;
    handles.ions_y = loaded_data.ions_y;
    handles.ions_tof = loaded_data.ions_tof;
    handles.shutterstatusRaw = loaded_data.shutterstatusRaw;
    handles.delaystatusRaw = loaded_data.delaystatusRaw;
    handles.polarizationstatusRaw = loaded_data.polarizationstatusRaw;
    handles.paramstatusRaw = loaded_data.paramstatusRaw;
    handles.intensitystatusRaw = loaded_data.intensitystatusRaw;
    handles.delayInfo = loaded_data.delayInfo;
    handles.polarInfo = loaded_data.polarInfo;
    handles.paramInfo = loaded_data.paramInfo;
    
    delayInfo = handles.delayInfo;
    x = cell(size((delayInfo(1,:))',1)+1,1);
    for nn = 1:size((delayInfo(1,:))',1)
        x(nn+1) = {num2str(delayInfo(1,nn))};
    end
    x(1) = {'all'};

    set(handles.delayChoice, 'string', x);
    
    polarInfo = handles.polarInfo;
    x = cell(size((polarInfo(1,:))',1)+1,1);
    for nn = 1:size((polarInfo(1,:))',1)
        x(nn+1) = {num2str(polarInfo(1,nn))};
    end
    x(1) = {'all'};

    set(handles.polarChoice, 'string', x);
    
    polarInfo = handles.paramInfo;
    x = cell(size(polarInfo(2,:),2)+1,1);
    x(2:end) = table2cell(polarInfo(2,:));
    x(1) = {'all'};
    x = lower(x);
    
    set(handles.paramChoice, 'string', x);
    
    clear loaded_data

function saveLoadedData(handles, delayInfo, polarInfo)
    loaded_data = struct('datafile', handles.datafile, 'path', handles.path,...
                         'numHitsRaw', handles.numHitsRaw, ... %  'HitNoRaw', handles.HitNoRaw,...
                         'ions_x', handles.ions_x, 'ions_y', handles.ions_y,...
                         'ions_tof', handles.ions_tof, 'closedshutterRaw', handles.shutterstatusRaw,...
                         'chosenDelayRaw', handles.delaystatusRaw, 'chosenPolarizationRaw',...
                         handles.polarizationstatusRaw, 'chosenParamRaw', handles.paramstatusRaw,...
                         'delayInfo', delayInfo, 'polarInfo', polarInfo);

    if ~exist([handles.path '\analysis'], 'dir')
        mkdir([handles.path '\analysis']);
        pause(1)
    end
    
    timee = clock;
    filename = [handles.datafile, '-loaded_data-', date, '-', num2str(timee(4)), '-',...
        num2str(timee(5)), '-', num2str(floor(timee(6))), '.mat'];
    filename = fullfile([handles.path '\analysis\'], filename);
    
    save(filename, 'loaded_data', '-v7.3');

    
% --- Executes on button press in findTof.
%finds the tof for a particle with some mass, charge and eV
function findTof_Callback(hObject, eventdata, handles)
% hObject    handle to findTof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set the error bar to no errors
set(handles.err, 'string', 'no errors');
set(handles.findTof, 'string', 'calculating');
drawnow

%get the necessary data from UI
commonParams = get(handles.commonParams, 'data');
calibParams = get(handles.calibParams, 'data');

%get the necessary data from UI vectors
V1 = commonParams(4);
VM = commonParams(5);
ss = commonParams(6);
mass = calibParams(1);
charge = calibParams(2);
eV = calibParams(3);

%look for potential errors
if charge == 0
    set(handles.err, 'string', 'entered charge is zero');
    error('entered charge is zero')
elseif mass == 0
    set(handles.err, 'string', 'entered mass is zero');
    error('entered mass is zero')
end    

%enter the default values
if V1 == 1, V1 = 2000; end
if VM == 1, VM = -2250; end
if ss == 1, ss = 40.5; end

%retreive flym data and display tof from simion
if eV > 0
    evalc('Sim = Flym_Sim(charge, mass, eV, 0, 0, ss, V1, VM);');
else
    evalc('Sim = Flym_Sim(charge, mass, abs(eV), 0, 180, ss, V1, VM);');
end

tof = Sim(2, 2)*10^3;
set(handles.tofNs, 'string', [num2str(tof, 5) ' (ns)']);

%save any values saved to handles
guidata(hObject, handles);

set(handles.findTof, 'string', 'find tof');


% --- Executes on button press in tofHist.
%creates histogram of tof or mass
function tofHist_Callback(hObject, eventdata, handles)
% hObject    handle to tofHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%reset the error bar and the button string
set(handles.err, 'string', 'no errors');
set(handles.tofHist, 'string', 'plotting');
drawnow

%get the necessary data from UI
tofNumBins = str2double(get(handles.tofNumBins, 'string'));
t0 = get(handles.commonParams, 'data');
calibPoints = get(handles.calibPoints, 'data');

%get the data
tof = handles.ions_tof;
rX = handles.ions_x;
rY = handles.ions_y;

%get the necessary data from UI vectors
calibMass = calibPoints(1, :);
calibCharge = calibPoints(2, :);
calibTofMin = calibPoints(3, :);
calibTofMax = calibPoints(4, :);

%check the tof and position are loaded
if ~isfield(handles, 'ions_tof')||~isfield(handles, 'ions_x')||~isfield(handles, 'ions_y')
    set(handles.err, 'string', 'data is not loaded');
    set(handles.tofHist, 'string', 'histograms');
    error('data is not loaded')
elseif tofNumBins == 0||isnan(tofNumBins)
    set(handles.err, 'string', 'number of bins is zero');
    set(handles.tofHist, 'string', 'histograms');
    error('number of bins is zero')
elseif prod(isnan(calibTofMin))||prod(isnan(calibTofMax))
    set(handles.err, 'string', 'Tof Min or Tof Max is not a number');
    set(handles.tofHist, 'string', 'histograms');
    error('Tof Min or Tof Max is not a number');
end

%set conditions for shutterChoice and paramChoice settings
cond = (tof(:, 1) > 50)&(tof(:, 2) > 50);
cond = ...
    ApplyExperimentType(cond, get(handles.shutterChoice, 'value'), ...
                        get(handles.intensityChoice, 'value'), get(handles.paramChoice, 'value'), ...
                        get(handles.polarChoice, 'value'), get(handles.delayChoice, 'value'),...
                        handles.shutterstatusRaw,   handles.intensitystatusRaw,...
                        handles.paramstatusRaw, handles.polarizationstatusRaw,...
                        handles.delaystatusRaw, handles.polarInfo, handles.delayInfo,...
                        handles.paramInfo);

%apply the conditions
tof = tof(cond, :);
rX = rX(cond, :);
rY = rY(cond, :);

plotPIPICO(tof, t0(1), 'tof 1st (ns)', 'tof 2st (ns)')

%plot histograms of tof, x position, and y position
% Exclude entries that are not events (assuming no ion has tof exactly 0)
plotCalibrationHist(rX(tof ~= 0),  t0(1), tofNumBins, 'x (mm)', true);
plotCalibrationHist(rY(tof ~= 0),  t0(2), tofNumBins, 'y (mm)', true);
[bins, values] = plotCalibrationHist(tof(tof ~= 0), t0(3), tofNumBins, 'tof (ns)', true);
hold on

%plot figures with the histogram tof fit
if prod(calibTofMin)*prod(calibTofMax) > 0
    [calibPoints(6,:)] = fitTOF(bins, values, calibTofMin, calibTofMax, calibMass, calibCharge);
    set(handles.calibPoints, 'data', calibPoints)
end
hold off

if get(handles.saveCalib, 'value')
    saveCalibration(handles, calibPoints)
end

%reset the button string
set(handles.tofHist, 'string', 'histograms');

%save any values saved to handles
guidata(hObject, handles);


function plotBrokenDataTof ()
    [loadFile, loadPath] = uigetfile(handles.path);
    loadFile = loadFile(1:end-5);
    allfiles = dir([loadPath, '\*.mat']);
    numfiles = length(allfiles);
    figure()
    for nn = 1:numfiles
        load([loadPath, loadFile, num2str(nn), '.mat']);
        tof = loaded_data.ions_tof;
        cond = (tof(:, 1) > 50)&(tof(:, 2) > 50);
        cond = ...
            ApplyExperimentType(cond, get(handles.shutterChoice, 'value'), ...
                                get(handles.intensityChoice, 'value'), get(handles.paramChoice, 'value'), ...
                                get(handles.polarChoice, 'value'), get(handles.delayChoice, 'value'),...
                                loaded_data.shutterstatusRaw,   loaded_data.intensitystatusRaw,...
                                loaded_data.paramstatusRaw, loaded_data.polarizationstatusRaw,...
                                loaded_data.delaystatusRaw, loaded_data.polarInfo, loaded_data.delayInfo,...
                                loaded_data.paramInfo);
        sum(cond)
        tof = tof(cond, :);
        
        plotCalibrationHist(tof(tof ~= 0), t0(3), tofNumBins, 'tof (ns)', false, ...
                            'Color', [1 - (nn-1)/numfiles, 0, (nn-1)/numfiles]);
        hold on 
    end
    hold off


function plotPIPICO(tof, t0, label1, label2)
if (t0 ~= 1)&&(t0 ~= 0)&& ~(isnan(t0))
    figure()
    plot(tof(:, 1)-t0, tof(:, 2)-t0, '.')
else
    figure()
    plot(tof(:, 1), tof(:, 2), '.')
end
xlabel(label1)
ylabel(label2)


function [bins, values] = plotCalibrationHist(xall, x0, num_bins, labelx, ...
                                              newfigure, varargin)
if (x0 ~= 1)&&~(isnan(x0))
    xall = xall-x0;
    [values, bins] = hist(xall, num_bins);
else
    [values, bins] = hist(xall, num_bins);
end
if (newfigure)
    figure()
end

plot(bins, values, varargin{:})
xlabel(labelx)
title(['number of hits ' num2str(numel(xall))])


function [tof_fits] = fitTOF(bins, values, calibTofMin, calibTofMax, ...
                                calibMass, calibCharge)
tof_peaks = [0, 0, 0];
colorArray = [1, .5, .5; .5, .75, .5; .75, 0, .75];

for ii = 1:3
    fit_range = (calibTofMin(ii) < bins)&(bins < calibTofMax(ii));
    try
        fit_result = fit(bins(fit_range).', values(fit_range).', 'gauss1');
        tof_peaks(ii) = fit_result.b1;
        p = plot(bins(fit_range), fit_result(bins(fit_range)));
        set(p, 'color', colorArray(ii, :))
    catch
        tof_peaks(ii) = NaN;
    end
end

legend('data', ['mass ' num2str(calibMass(1)) ' charge '...
    num2str(calibCharge(1)) ' tof peak ' num2str(tof_peaks(1))],...
    ['mass ' num2str(calibMass(2)) ' charge '...
    num2str(calibCharge(2)) ' tof peak ' num2str(tof_peaks(2))],...
    ['mass ' num2str(calibMass(3)) ' charge '...
    num2str(calibCharge(3)) ' tof peak ' num2str(tof_peaks(3))])

tof_fits = tof_peaks;



function saveCalibration(handles, calibPoints)
    points_of_interest = struct('calibPoints', calibPoints, 'filename', handles.datafile);

    if ~exist([handles.path '\analysis'], 'dir')
        mkdir([handles.path '\analysis']);
        pause(1)
    end
    
    timee = clock;
    filename = [handles.datafile, '-points_of_interest-', date, '-', num2str(timee(4)), '-',...
        num2str(timee(5)), '-', num2str(floor(timee(6))), '.mat'];
    filename = fullfile([handles.path '\analysis\'], filename);
    
    save(filename, 'points_of_interest', '-v7.3');




function tweekParams_Callback(hObject, eventdata, handles)
% hObject    handle to tweekParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%reset the error bar and the button string
set(handles.err, 'string', 'no errors');
set(handles.tweekParams, 'string', 'plotting');
pause(1)

%get the necessary data from UI
commonParams = get(handles.commonParams, 'data');
saveCalib = get(handles.saveCalib, 'value');
calibPoints = get(handles.calibPoints, 'data');
shutterChoice = get(handles.shutterChoice, 'value');
intensityChoice = get(handles.intensityChoice, 'value');
polarChoice = get(handles.polarChoice, 'value');
paramChoice = get(handles.paramChoice, 'value');
delayChoice = get(handles.delayChoice, 'value');

%get the data
tof = handles.ions_tof;
rX = handles.ions_x;
rY = handles.ions_y;

%get the necessary data from UI vectors
x0 = commonParams(1);
y0 = commonParams(2);
t0 = commonParams(3);
V1 = commonParams(4);
VM = commonParams(5);
ss = commonParams(6);
calibMass = calibPoints(1, :);
calibCharge = calibPoints(2, :);
calibTofMin = calibPoints(3, :);
calibTofMax = calibPoints(4, :);
calibTof = calibPoints(5, :);

%check for potential errors
if ~isfield(handles, 'ions_tof')||~isfield(handles, 'ions_x')||~isfield(handles, 'ions_y')
    set(handles.err, 'string', 'data is not loaded');
    set(handles.tweekParams, 'string', 'calibration');
    error('data is not loaded')
elseif ~(sum(calibMass) > 2)||~(sum(calibCharge) > 2)||~(sum(calibTof) > 2)
    set(handles.err, 'string', '3 points must be chosen to do the calibraion');
    set(handles.tweekParams, 'string', 'calibration');
    error('3 points must be chosen to do the calibraion');
elseif (calibTofMax(1)*calibTofMin(1) < 1)&&(x0 == 1)&&(y0 == 1)
    set(handles.err, 'string', 'enter a min max tof for your first species to find x0 and y0');
    set(handles.tweekParams, 'string', 'calibration');
    error('enter a min max tof for your first species to find x0 and y0');
end

%set conditions
cond = (tof(:, 1) > 50)&(tof(:, 2) > 50)&(tof(:, 3) > 50)&(tof(:, 4) > 50);
%'main line 742, set the previous line back on'
%cond = (tof(:, 1) > 50);

cond = ApplyExperimentType(cond, shutterChoice, intensityChoice,...
    paramChoice, polarChoice, delayChoice, handles.shutterstatusRaw, handles.intensitystatusRaw, handles.paramstatusRaw,...
    handles.polarizationstatusRaw, handles.delaystatusRaw, handles.polarInfo, handles.delayInfo,...
    handles.paramInfo);

%load the file and then get the positions and tof
tof = tof(cond, :);
rX = rX(cond, :);
rY = rY(cond, :);
tof = tof(:);
rX = rX(:);
rY = rY(:);

%take care of default values
if V1 == 1, V1 = 2000; end
if VM == 1, VM = -2250; end

%find the x0 and y0
if x0 == 1
    x0 = mean(rX(tof > calibTofMin(1) & tof < calibTofMax(1)));
end
if y0 == 1
    y0 = mean(rY(tof > calibTofMin(1) & tof < calibTofMax(1)));
end

tofNumBins = str2double(get(handles.tofNumBins, 'string'));
figure();
hist(rX - x0, tofNumBins);
hold on;
hist(rX(tof > calibTofMin(1) & tof < calibTofMax(1)) - x0, tofNumBins);
xlabel('x (mm)')
title('Distribution of all and calibration species')


figure();
hist(rY - y0, tofNumBins);
hold on;
hist(rY(tof > calibTofMin(1) & tof < calibTofMax(1)) - y0, tofNumBins);
xlabel('y (mm)')
title('Distribution of all and calibration species')

%check for some errors
if V1 < 500 || VM > -500
    set(handles.err, 'string',...
        'voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500');
    set(handles.tweekParams, 'string', 'calibration');
    error('voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500')
end

%now do the calibration if ss is left as default
if ss == 1
    %make possible initial positions
    posTOF = linspace(36.5,43,100);  %x, along TOF axis, 12.5 is center
    N_pos = length(posTOF);

    %retreive fly'm data from simion
    Sim = Flym_Sim(calibCharge, calibMass, 0, 0, 0, posTOF, V1, VM);
    t_Sim = Sim(2:2:end, 2)*10^3;
    
    mm = 0;
    
    %if t0 was default, calibrate this too
    if t0 == 1
        %make possible t0 values and initialize necessary arrays
        t0_trying = linspace(0, 200, 500);
        Diffs = zeros(N_pos*length(t0_trying),1);
        II = Diffs;
        JJ = Diffs;
        
        %find difference between calculated tof and actual to minimize later
        for jj = 1:length(t0_trying)
            for ii = 0:N_pos-1
                mm = mm+1;

                Diffs(mm) = abs(t_Sim(3*ii+1)-calibTof(1)+t0_trying(jj))+...
                    abs(t_Sim(3*ii+2)-calibTof(2)+t0_trying(jj))+...
                    abs(t_Sim(3*ii+3)-calibTof(3)+t0_trying(jj));

                II(mm) = ii;
                JJ(mm) = jj;
            end
            [~, I_min] = min(Diffs);
        end
        
        %retrieve the calibrated t0
        calibrated_t0 = t0_trying(JJ(I_min));

    else
        %calibrate t0 and ss
        %start by initializing necessary arrays
        Diffs = zeros(N_pos,1);
        II = Diffs;
        
        %find difference between calculated tof and actual to minimize later
        for ii = 0:N_pos-1
            mm = mm+1;

            Diffs(mm) = abs(t_Sim(3*ii+1)-calibTof(1)+t0)+abs(t_Sim(3*ii+2)-calibTof(2)+t0)+...
                abs(t_Sim(3*ii+3)-calibTof(3)+t0);

            II(mm) = ii;
        end
        [~, I_min] = min(Diffs);
        
        %t0 was given by user input
        calibrated_t0 = t0;
    end

    %retrieve the calibrated ss
    calibrated_s = posTOF(II(I_min)+1);
    
    %plot user input tof
    C_calc = [0.7 0.3 0.3];
    C_sim = [0.3 0.3 0.7];
    figure()
    plot(calibMass./calibCharge, calibTof-calibrated_t0, 'Marker', 'x', 'LineStyle', 'none', 'Color', C_sim)
    
    %plot calculated tof after calibration if it was required
    if t0 == 1||ss == 1
        hold on
        plot(calibMass./calibCharge, [t_Sim(3*(II(I_min))+1), t_Sim(3*(II(I_min))+2),...
            t_Sim(3*(II(I_min))+3)], 'Marker', 'o', 'LineStyle', 'none', 'Color', C_calc)
        xlabel('mass/charge')
        ylabel('t (ns)')
        hold off
    end
   
%now do the calibration if ss is a user input
elseif ~(ss == 1)
    t_Sim = zeros(3, 1);
    
    %calculate the tof and output the value to UI
    for ii = 1:3
        %retreive fly'm data from Simion
        Sim = Flym_Sim(calibCharge(ii), calibMass(ii), 0, 0, 0, ss, V1, VM);
        t_Sim(ii) = Sim(2:2:end, 2)*10^3;
    end
    
    %make array of possible t0 and calculate difference between actual tof
    %and calculated tof to minimize
    t0_trying = linspace(0, 200, 500);
    Diffs = abs(t_Sim(1)-calibTof(1)+t0_trying)+abs(t_Sim(2)-calibTof(2)+t0_trying)+...
        abs(t_Sim(3)-calibTof(3)+t0_trying);
    [~, I_min] = min(Diffs);
    
    %set s as the user input value and generate calibrated t0 value
    calibrated_s = ss;
    if t0 == 1
        calibrated_t0 = t0_trying(I_min);
    else
        calibrated_t0 = t0;
    end
    
    %plot best one
    C_calc = [0.7 0.3 0.3];
    C_sim = [0.3 0.3 0.7];
    figure()
    plot(calibMass./calibCharge, calibTof-calibrated_t0, 'Marker', 'x', 'LineStyle', 'none', 'Color', C_sim)
    hold on
    plot(calibMass./calibCharge, [t_Sim(1), t_Sim(2), t_Sim(3)],...
        'Marker', 'o', 'LineStyle', 'none', 'Color', C_calc)
    xlabel('mass/charge')
    ylabel('t (ns)')
    hold off
end

%if no calibration was necessary, use user input values
if ~(t0 == 1)&&~(ss == 1)
    calibrated_t0 = t0;
    calibrated_s = ss;
end

% Export starting position
set(handles.s_fit, 'string', ['s = ' num2str(calibrated_s) ' mm']);
set(handles.t0_fit, 'string', ['t0 = ' num2str(calibrated_t0) ' ns']);
set(handles.x0_fit, 'string', ['x0 = ' num2str(x0) ' mm']);
set(handles.y0_fit, 'string', ['y0 = ' num2str(y0) ' mm']);

%save outputs
if saveCalib
    calib = struct('x0', x0, 'y0', y0, 't0', calibrated_t0, 'V1', V1,...
        'VM', VM, 's', calibrated_s, 'calibration', calibPoints, 'filename', handles.datafile);
    
    if ~exist([handles.path '\analysis'], 'dir')
        mkdir([handles.path '\analysis']);
        pause(1)
    end
    
    timee = clock;
    filename = [handles.datafile, '-calibration-', date, '-', num2str(timee(4)), '-',...
        num2str(timee(5)), '-', num2str(floor(timee(6))), '.mat'];
    filename = fullfile([handles.path '\analysis\'], filename);
    
    save(filename, 'calib');

end

%reset the button string
set(handles.tweekParams, 'string', 'calibration');

%save any values saved to handles
guidata(hObject, handles);



% --- Executes on button press in prepare.
%gets everything together from the GUI so that the prepare function can be
%called
function prepare_Callback(hObject, eventdata, handles)
% hObject    handle to prepare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%reset the error bar and the button string
set(handles.err, 'string', 'no errors');
set(handles.prepare, 'string', 'preparing');
pause(1)

%get the necessary data from UI
commonParams = get(handles.commonParams, 'data');
prepareParams = get(handles.prepareParams, 'data');
saveData = get(handles.savePrepare, 'value');
includePrepared = get(handles.includePrepared, 'value');
loadCalib = get(handles.loadCalib, 'value');
useBrokeData = get(handles.useBrokeData, 'value');
calibPoints = get(handles.calibPoints, 'data');
autoCalib = get(handles.autoCalib, 'value');

%set length of ev array and theta array
EVlength = 25;
Thetalength = 80;

%get the necessary data from UI vectors
x0 = commonParams(1);
y0 = commonParams(2);
t0 = commonParams(3);
V1 = commonParams(4);
VM = commonParams(5);
ss = commonParams(6);
mass = prepareParams(1, :);
charge = prepareParams(2, :);
maxEV = prepareParams(3, :);
mass = mass(mass > 0);
charge = charge(charge > 0);
maxEV = maxEV(maxEV > 0);

%check for potential errors
if length(mass) ~= length(charge)
    set(handles.err, 'string', 'number of masses does not equal number of charges');
    set(handles.prepare, 'string', 'prepare');
    error('number of masses does not equal number of charges')
elseif isempty(mass)
    set(handles.err, 'string', 'masses are empty');
    set(handles.prepare, 'string', 'prepare');
    error('masses are empty')
elseif isempty(charge)
    set(handles.err, 'string', 'charge are empty');
    set(handles.prepare, 'string', 'prepare');
    error('charge are empty')
end

if loadCalib

    set(handles.err, 'string', 'please select your calibration file');
    
    if exist('handles', 'var') && isfield(handles, 'path')
        [calibFile, calibPath] = uigetfile(handles.path);
    else
        [calibFile, calibPath] = uigetfile();
    end
    
    set(handles.err, 'string', 'no errors');
    
    load([calibPath, calibFile])
    
    x0 = calib.x0;
    y0 = calib.y0;
    t0 = calib.t0;
    V1 = calib.V1;
    VM = calib.VM;
    ss = calib.s;
    
    set(handles.commonParams, 'data', [x0, y0, t0, V1, VM, ss]);
    drawnow
    
end

if includePrepared

    if useBrokeData
        set(handles.err, 'string', 'you cannot include previously prepared data with broken up data')
        error('you cannot include previously prepared data with broken up data')
    end
    
    set(handles.err, 'string', 'please select your prepared file');
    
    if exist('handles', 'var') && isfield(handles, 'path')
        [calibFile, calibPath] = uigetfile(handles.path);
    else
        [calibFile, calibPath] = uigetfile();
    end
    
    set(handles.err, 'string', 'no errors');
    
    load([calibPath, calibFile])
    
    notincluded = [];
    notindex = 0;
    preparedIndex = [];
    
    for ii = 1:length(mass)

        index = find(mass(ii) == prepared.mass & charge(ii) == prepared.charge & maxEV(ii) == prepared.maxEV);
        
        if isempty(index)
            notindex = notindex+1;
            notincluded(notindex,:) = [mass(ii), charge(ii), maxEV(ii)];
        else
            preparedIndex(ii-notindex) = find(mass(ii) == prepared.mass & charge(ii) == prepared.charge & maxEV(ii) == prepared.maxEV);
        end
    end
    
    if notindex == 0
        
        handles.EV = prepared.EV(:, preparedIndex);
        handles.momZ = prepared.momZ(:, preparedIndex);
        handles.momX = prepared.momX(:, preparedIndex);
        handles.momY = prepared.momY(:, preparedIndex);
        handles.hitNo = prepared.hitNo;
        handles.shotNo = prepared.shotNo;
        handles.numHits_processed = prepared.numHits_processed;
        handles.ions_tof_processed = prepared.ions_tof_processed;
        handles.ions_x_processed = prepared.ions_x_processed;
        handles.ions_y_processed = prepared.ions_y_processed;
        handles.shutterStatus = prepared.shutterStatus;
        handles.intensityStatus = prepared.intensityStatus;
        handles.polarizationStatus = prepared.polarizationStatus;
        handles.paramStatus = prepared.paramStatus;
        handles.delayStatus = prepared.delayStatus;
        
    else
        
        if ~isfield(handles, 'ions_tof')||~isfield(handles, 'ions_x')||~isfield(handles, 'ions_y')
            set(handles.err, 'string', 'data is not loaded');
            set(handles.prepare, 'string', 'prepare');
            error('data is not loaded')
        end
    
        %take care of default values
        if x0 == 1, x0 = 0; end
        if y0 == 1, y0 = 0; end
        if t0 == 1, t0 = 0; end
        if V1 == 1, V1 = 2000; end
        if VM == 1, VM = -2250; end
        if ss == 1, ss = 40.5; end

        if V1 < 500 || VM > -500 && (V1 ~= 1) && (VM ~= 1)
            set(handles.err, 'string',...
                'voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500');
            set(handles.prepare, 'string', 'prepare');
            error('voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500')
        end
       
        [eVArray, thetaArray, tof_Sim, r_Sim] = makeSimArrays(V1, VM, ss, mass, charge, maxEV, EVlength, Thetalength);
        
        [handles.EV, handles.momZ, handles.momX, handles.momY, handles.hitNo, handles.shotNo,...
        handles.numHits_processed, handles.ions_tof_processed, handles.ions_x_processed, handles.ions_y_processed,...
        handles.shutterStatus, handles.intensityStatus, handles.polarizationStatus, handles.paramStatus, handles.delayStatus] =...
        prepare(V1, VM, ss, t0, x0, y0, notincluded(:, 1), notincluded(:, 2), notincluded(:, 3), handles.ions_tof,...
        handles.ions_x, handles.ions_y, handles.numHitsRaw, eVArray, thetaArray, tof_Sim, r_Sim,...
        handles.shutterstatusRaw, handles.intensityStatusRaw, handles.polarizationStatusRaw,...
        handles.paramStatusRaw, handles.delayStatusRaw, 0);
        
        handles.EV = [handles.EV, prepared.EV(:, preparedIndex)];
        handles.momZ = [handles.momZ, prepared.momZ(:, preparedIndex)];
        handles.momX = [handles.momX, prepared.momX(:, preparedIndex)];
        handles.momY = [handles.momY, prepared.momY(:, preparedIndex)];
        
        mass = [notincluded(:, 1), prepared.mass(preparedIndex)];
        charge = [notincluded(:, 2), prepared.charge(preparedIndex)];
        maxEV = [notincluded(:, 3), prepared.maxEV(preparedIndex)];

        prepareParams = zeros(size(handles.prepareParams));
        prepareParams([1,2,3], linspace(1,length(mass),length(mass))) = [mass;charge;maxEV];
        set(handles.prepareParams, 'data', prepareParams)
        
    end
    
    clear prepared
    
elseif useBrokeData

    set(handles.err, 'string', 'please select one loaded file');
    
    if exist('handles', 'var') && isfield(handles, 'path')
        [loadFile, loadPath] = uigetfile(handles.path);
    else
        [loadFile, loadPath] = uigetfile();
        handles.path = loadPath;
    end
    
    set(handles.err, 'string', 'no errors');
    
    loadFile = loadFile(1:end-5);
    handles.datafile = loadFile;
    
    numfiles = dir([loadPath, '\*.mat']);
    numfiles = length(numfiles);

    %take care of default values
    if x0 == 1, x0 = 0; end
    if y0 == 1, y0 = 0; end
    if t0 == 1, t0 = 0; end
    if V1 == 1, V1 = 2000; end
    if VM == 1, VM = -2250; end
    if ss == 1, ss = 40.5; end

    if V1 < 500 || VM > -500 && (V1 ~= 1) && (VM ~= 1)
        set(handles.err, 'string',...
            'voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500');
        set(handles.prepare, 'string', 'prepare');
        error('voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500')
    end
    
    %call the prepare function that will generate the momentum matrix that is requried

    handles.EV = [];
    handles.momZ = [];
    handles.momX = [];
    handles.momY = [];
    handles.hitNo = [];
    handles.shotNo = [];
    handles.numHits_processed = [];
    handles.ions_tof_processed = [];
    handles.ions_x_processed = [];
    handles.ions_y_processed = [];
    handles.shutterStatus = [];
    handles.intensityStatus = [];
    handles.polarizationStatus = [];
    handles.paramStatus = [];
    handles.delayStatus = [];
    
    
    [eVArray, thetaArray, tof_Sim, r_Sim] = makeSimArrays(V1, VM, ss, mass, charge, maxEV, EVlength, Thetalength);
    
    ss_old = str2double(erase(erase(get(handles.s_fit, 'string'), 's = '), ' mm'));
    t0_old = str2double(erase(erase(get(handles.t0_fit, 'string'), 't0 = '), ' ns'));
    x0_old = str2double(erase(erase(get(handles.x0_fit, 'string'), 'x0 = '), ' mm'));
    y0_old = str2double(erase(erase(get(handles.y0_fit, 'string'), 'y0 = '), ' mm'));
    
    tic    
    
    for nn = 1:numfiles
 
        disp(['working on file number ' num2str(nn) ' of ' num2str(numfiles)])
        
        load([loadPath, loadFile, num2str(nn), '.mat']);
        
        if autoCalib
            
            mainfigurehandle = gcf;
            
            handles.ions_x = loaded_data.ions_x;
            handles.ions_y = loaded_data.ions_y;
            handles.ions_tof = loaded_data.ions_tof;
            
            calibPoints(6,:) = [0,0,0];
            set(handles.calibPoints, 'data', calibPoints)
            commonParams(3) = 0;
            set(handles.commonParams, 'data', commonParams)
            tofHist_Callback(hObject, eventdata, handles)
            calibPoints = get(handles.calibPoints, 'data');
            calibPoints(5,:) = calibPoints(6,:);
            set(handles.calibPoints, 'data', calibPoints)
            
            commonParams([1,2,3,6]) = 1;
            set(handles.commonParams, 'data', commonParams)
            
            try
                tweekParams_Callback(hObject, eventdata, handles)

                ss = str2double(erase(erase(get(handles.s_fit, 'string'), 's = '), ' mm'));
                t0 = str2double(erase(erase(get(handles.t0_fit, 'string'), 't0 = '), ' ns'));
                x0 = str2double(erase(erase(get(handles.x0_fit, 'string'), 'x0 = '), ' mm'));
                y0 = str2double(erase(erase(get(handles.y0_fit, 'string'), 'y0 = '), ' mm'));
                
                ss_old = ss;
                t0_old = t0;
                
                if isnan(x0)
                    x0 = x0_old;
                end
                
                if isnan(y0)
                    y0 = y0_old;
                end
                
                x0_old = x0;
                y0_old = y0;
            
            catch
                
                ss = ss_old;
                t0 = t0_old;
                
                if isnan(x0)
                    x0 = x0_old;
                end
                
                if isnan(y0)
                    y0 = y0_old;
                end
                
                x0_old = x0;
                y0_old = y0;

            end
            

            
            [eVArray, thetaArray, tof_Sim, r_Sim] = makeSimArrays(V1, VM, ss, mass, charge, maxEV, EVlength, Thetalength);
            
            set(mainfigurehandle, 'HandleVisibility', 'off');
            close all
            set(mainfigurehandle, 'HandleVisibility', 'on');
            
        end
        
        [EV, momZ, momX, momY, hitNo, shotNo, numHits_processed, ions_tof_processed,...
            ions_x_processed, ions_y_processed, shutterStatus, intensityStatus,...
            polarizationStatus, paramStatus, delayStatus] =...
            prepare(V1, VM, ss, t0, x0, y0, mass, charge, maxEV, loaded_data.ions_tof, loaded_data.ions_x,...
            loaded_data.ions_y, loaded_data.numHitsRaw, eVArray, thetaArray, tof_Sim, r_Sim,...
            loaded_data.shutterstatusRaw, loaded_data.intensitystatusRaw, loaded_data.polarizationstatusRaw,...
            loaded_data.paramstatusRaw, loaded_data.delaystatusRaw, loaded_data.shotsStartZero);
        
        handles.EV = [handles.EV; EV];
        handles.momZ = [handles.momZ; momZ];
        handles.momX = [handles.momX; momX];
        handles.momY = [handles.momY; momY];
        handles.hitNo = [handles.hitNo; hitNo];
        handles.shotNo = [handles.shotNo; shotNo];
        handles.numHits_processed = [handles.numHits_processed; numHits_processed];
        handles.ions_tof_processed = [handles.ions_tof_processed; ions_tof_processed];
        handles.ions_x_processed = [handles.ions_x_processed; ions_x_processed];
        handles.ions_y_processed = [handles.ions_y_processed; ions_y_processed];
        handles.shutterStatus = [handles.shutterStatus; shutterStatus];
        handles.intensityStatus = [handles.intensityStatus; intensityStatus];
        handles.polarizationStatus = [handles.polarizationStatus; polarizationStatus];
        handles.paramStatus = [handles.paramStatus; paramStatus];
        handles.delayStatus = [handles.delayStatus; delayStatus];
    
        estimatedTime = toc/nn*numfiles - toc;
        
        set(handles.err, 'string', ['estimated completed load time: '...
        datestr(datetime('now')+seconds(estimatedTime))]);
        drawnow
        
    end
    
    handles.shotsStartZero = loaded_data.shotsStartZero;
    handles.polarInfo = loaded_data.polarInfo;
    handles.delayInfo = loaded_data.delayInfo;
    
else
    
    if ~isfield(handles, 'ions_tof')||~isfield(handles, 'ions_x')||~isfield(handles, 'ions_y')
        set(handles.err, 'string', 'data is not loaded');
        set(handles.prepare, 'string', 'prepare');
        error('data is not loaded')
    end

    %take care of default values
    if x0 == 1, x0 = 0; end
    if y0 == 1, y0 = 0; end
    if t0 == 1, t0 = 0; end
    if V1 == 1, V1 = 2000; end
    if VM == 1, VM = -2250; end
    if ss == 1, ss = 40.5; end

    if V1 < 500 || VM > -500 && (V1 ~= 1) && (VM ~= 1)
        set(handles.err, 'string',...
            'voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500');
        set(handles.prepare, 'string', 'prepare');
        error('voltage values are unreasonable, V1 should be greater than 500 and VM should be less than -500')
    end

    %call the prepare function that will generate the momentum matrix that is requried
    
    [eVArray, thetaArray, tof_Sim, r_Sim] = makeSimArrays(V1, VM, ss, mass, charge, maxEV, EVlength, Thetalength);
    
    [handles.EV, handles.momZ, handles.momX, handles.momY, handles.hitNo, handles.shotNo,...
        handles.numHits_processed, handles.ions_tof_processed, handles.ions_x_processed, handles.ions_y_processed,...
        handles.shutterStatus, handles.intensityStatus, handles.polarizationStatus,...
        handles.paramStatus, handles.delayStatus] =...
        prepare(V1, VM, ss, t0, x0, y0, mass, charge, maxEV, handles.ions_tof, handles.ions_x,...
        handles.ions_y, handles.numHitsRaw, eVArray, thetaArray, tof_Sim, r_Sim,...
        handles.shutterstatusRaw, handles.intensitystatusRaw, handles.polarizationstatusRaw,...
        handles.paramStatusRaw, handles.delayStatusRaw, 0);

end

if saveData
    %build the output file

    if ~exist([handles.path '\analysis'], 'dir')
        mkdir([handles.path '\analysis']);
        pause(1)
    end
    
    timee = clock;
    filename = [handles.datafile, '-prepared-', date, '-', num2str(timee(4)), '-',...
        num2str(timee(5)), '-', num2str(floor(timee(6))), '.mat'];
    filename = fullfile([handles.path '\analysis\'], filename);
    
    prepared.EV = handles.EV;
    prepared.momZ = handles.momZ;
    prepared.momX = handles.momX;
    prepared.momY = handles.momY;
    prepared.hitNo = handles.hitNo;
    prepared.shotNo = handles.shotNo;
    prepared.numHits_processed = handles.numHits_processed;
    prepared.ions_tof_processed = handles.ions_tof_processed;
    prepared.ions_x_processed = handles.ions_x_processed;
    prepared.ions_y_processed = handles.ions_y_processed;
    prepared.shutterStatus = handles.shutterStatus;
    prepared.intensityStatus = handles.intensityStatus;
    prepared.polarizationStatus = handles.polarizationStatus;
    prepared.paramStatus = handles.paramStatus;
    prepared.delayStatus = handles.delayStatus;
    prepared.shotsStartZero = handles.shotsStartZero;
    prepared.polarInfo = handles.polarInfo;
    prepared.delayInfo = handles.delayInfo;
    prepared.mass = mass;
    prepared.charge = charge;
    prepared.maxEV = maxEV;
    prepared.datafile = handles.datafile;
    prepared.path = handles.path;

    save(filename, 'prepared', '-v7.3');
end

%reset the button string
set(handles.prepare, 'string', 'prepare');

%save any values saved to handles
guidata(hObject, handles);


% --- Executes on button press in cutAndPlot.
%gets everything together from the GUI so that the cutandplot function can
%be called
function cutAndPlot_Callback(hObject, eventdata, handles)
% hObject    handle to cutAndPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%reset the error bar and the button string
set(handles.err, 'string', 'no errors');
set(handles.cutAndPlot, 'string', 'plotting');
pause(1)

%get the necessary data from UI
numExtraHits = get(handles.numExtraHits, 'string');
coneV = get(handles.coneV, 'value');
conDeltat = get(handles.conDeltat, 'value');
prepareParams = get(handles.prepareParams, 'data');
cutAndPlotParams = get(handles.cutAndPlotParams, 'data');
deltat = get(handles.deltat, 'string');
deltax = get(handles.deltax, 'string');
deltay = get(handles.deltay, 'string');

plTofPlot = get(handles.plTofPlot, 'value');
plTofHist = get(handles.plTofHist, 'value');
plEnergyPartHist = get(handles.plEnergyPartHist, 'value');
plEnergyHist = get(handles.plEnergyHist, 'value');
plMomDir = get(handles.plMomDir, 'value');
numCutPlotBins = get(handles.numCutPlotBins, 'string');
saveData = get(handles.saveData, 'value');
singles = get(handles.singles, 'value');
doubles = get(handles.doubles, 'value');
deleteFiles = get(handles.deleteFiles, 'value');

%get the necessary data from UI vectors
mass = cutAndPlotParams(1, :);
charge = cutAndPlotParams(2, :);
maxeV = cutAndPlotParams(3, :);
incl = cutAndPlotParams(4, :);
mass = mass(mass ~= 0);
charge = charge(charge ~= 0);
maxeV = maxeV(maxeV ~= 0);
incl = incl(mass ~= 0);
massPrep = prepareParams(1, :);
chargePrep = prepareParams(2, :);
massPrep = massPrep(massPrep ~= 0);
chargePrep = chargePrep(chargePrep ~= 0);

%check for potential errors
if length(mass) ~= length(charge) || length(mass) ~= length(maxeV)...
        || length(maxeV) ~= length(charge)
    set(handles.err, 'string',...
        'the number of charges, mass, and or max eV are not the same');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('the number of charges, mass, and or max eV are not the same')
end

%convert UI strings to doubles after checking for errors
if isempty(numCutPlotBins) 
    numCutPlotBins = 0;
else
    numCutPlotBins = str2double(numCutPlotBins);
end

if isempty(numExtraHits) 
    numExtraHits = 0;
else
    numExtraHits = str2double(numExtraHits);
    
    if numExtraHits+length(mass) > 8
        set(handles.err, 'string', 'too many extra particles requested');
        set(handles.cutAndPlot, 'string', 'cut and plot');
        error('too many extra particles requested')
    end
end

if length(mass)+numExtraHits > 4
    set(handles.err, 'string',...
        'you are trying to find a coincidence with more than 4 particles');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('you are trying to find a coincidence with more than 4 particles')
end

if isempty(deltat)
    deltat = 0;
else
    deltat = str2double(deltat);
end

if isempty(deltax)
    deltax = 0;
else
    deltax = str2double(deltax);
end

if isempty(deltay)
    deltay = 0;
else
    deltay = str2double(deltay);
end

%check for more potential errors
if isnan(numExtraHits)
    set(handles.err, 'string', 'number of extra particles is not a number');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('number of extra particles is not a number')
elseif (plMomDir && length(mass)<2)
    set(handles.err, 'string', 'at least two particles are required to plot mom sums');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('at least two particles are required to plot mom sums')
elseif (plTofPlot && length(mass)<2)
    set(handles.err, 'string', 'at least two particles are required to make TOF plot');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('at least two particles are required to make TOF plot')
elseif isnan(numCutPlotBins)
    set(handles.err, 'string', 'number of bins for histograms is not a number');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('number of bins for histograms is not a number')
elseif isnan(deltat)
    set(handles.err, 'string', 'deltat is not a number');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('deltat is not a number')
elseif isnan(deltax)
    set(handles.err, 'string', 'deltax is not a number');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('deltax is not a number')
elseif isnan(deltay)
    set(handles.err, 'string', 'deltay is not a number');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('deltay is not a number')
elseif numCutPlotBins == 0 && (plTofHist || plEnergyPartHist || plEnergyHist)
    set(handles.err, 'string', 'number of bins for histograms is not set or zero');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('number of bins for histograms is not set or zero')
elseif sum(incl) < 1
    set(handles.err, 'string', 'no particles are included in the coincidence');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('no particles are included in the coincidence')
elseif (conDeltat && length(mass)<2)
    set(handles.err, 'string', 'at least two particles are required for momentum sums');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    error('at least two particles are required for momentum sums')
end

%setup the colms vector that will connect that masses input for plotting
%and the masses used in the preparation
colms = zeros(size(mass));
for i = 1:length(mass)
    var = find(massPrep == mass(i) & chargePrep == charge(i));

    if isempty(var)
        set(handles.err, 'string',...
            'species in the coincidence was not in the preparation');
        set(handles.cutAndPlot, 'string', 'cut and plot');
        error('species in the coincidence was not in the preparation')
    else
        colms(i) = var;
    end
end

if singles
    multiplicity = 1;
    Sextras = sum(~incl);
    Textras = Sextras+1;
elseif doubles
    multiplicity = 2;
    Sextras = sum(~incl);

    kk=1;
    for ii = 1:sum(~incl)
        for jj = ii:sum(~incl)
            II(kk) = ii;
            JJ(kk) = jj;
            kk=kk+1;
        end
    end
    
    Dextras = length(II);
    Textras = Sextras+Dextras+1;
  
else
    Textras = 1;
end

MASS = mass;
CHARGE = charge;
COLMS = colms;
INCL = incl;
MAXEV = maxeV;
EXTRAMASS = mass(INCL == 0);
EXTRACHARGE = charge(INCL == 0);
EXTRACOLMS = colms(INCL == 0);
EXTRAMAXEV = maxeV(INCL == 0);

includeExtras = false;
doDoubles = false;

for XX = 1:Textras

    if (singles || doubles) && ~includeExtras
        mass = MASS(INCL == 1);
        charge = CHARGE(INCL == 1);
        colms = COLMS(INCL == 1);
        incl = INCL(INCL == 1);
        maxeV = MAXEV(INCL == 1);
        includeExtras = true;

    elseif (singles && includeExtras) || (doubles && includeExtras) && ~doDoubles
        mass = [MASS(INCL == 1), EXTRAMASS(XX-1)];
        charge = [CHARGE(INCL == 1), EXTRACHARGE(XX-1)];
        colms = [COLMS(INCL == 1), EXTRACOLMS(XX-1)];
        incl = [INCL(INCL == 1), 0];
        maxeV = [MAXEV(INCL == 1), EXTRAMAXEV(XX-1)];
       
        Sextras = Sextras - 1;
        if Sextras == 0
            doDoubles = true;
            X0 = XX;
        end

    elseif doubles && includeExtras && doDoubles
        mass = [MASS(INCL == 1), EXTRAMASS(II(XX-X0)), EXTRAMASS(JJ(XX-X0))];
        charge = [CHARGE(INCL == 1), EXTRACHARGE(II(XX-X0)), EXTRACHARGE(JJ(XX-X0))];
        colms = [COLMS(INCL == 1), EXTRACOLMS(II(XX-X0)), EXTRACOLMS(JJ(XX-X0))];
        incl = [INCL(INCL == 1), 0, 0];
        maxeV = [MAXEV(INCL == 1), EXTRAMAXEV(II(XX-X0)), EXTRAMAXEV(JJ(XX-X0))];

    end
    
    %call the cutandplot function to make the necessary cuts and to plot
    %everything

    cond = true(size(handles.shutterStatus));
    
    shutterChoice = get(handles.shutterChoice, 'value');
    intensityChoice = get(handles.intensityChoice, 'value');
    paramChoice = get(handles.paramChoice, 'value');
    polarChoice = get(handles.polarChoice, 'value');
    delayChoice = get(handles.delayChoice, 'value');
    
    cond = ApplyExperimentType(cond, shutterChoice, intensityChoice,...
        paramChoice, polarChoice, delayChoice, handles.shutterStatus, handles.intensityStatus, handles.paramStatus,...
        handles.polarizationStatus, handles.delayStatus, handles.polarInfo, handles.delayInfo,...
        handles.paramInfo);
        
    [output] = cutandplot(mass, charge, colms, handles.momX(cond, :), handles.momY(cond, :),...
        handles.momZ(cond, :), handles.EV(cond, :), handles.ions_tof_processed(cond),...
        handles.hitNo(cond), handles.shotNo(cond), handles.numHits_processed(cond),...
        numExtraHits, coneV, maxeV, plTofPlot, plTofHist, plEnergyPartHist, plEnergyHist,...
        plMomDir, numCutPlotBins, deltat, deltax, deltay, conDeltat, incl);

    if saveData
        %build the output file
        
        if ~exist([handles.path '\analysis'], 'dir')
        mkdir([handles.path '\analysis']);
        pause(1)
        end

        timee = clock;
        filename = [handles.datafile, '-cut-', date, '-', num2str(timee(4)), '-',...
            num2str(timee(5)), '-', num2str(floor(timee(6))), '.mat'];
        filename = fullfile([handles.path '\analysis\'], filename);

        save(filename, 'output');
    end
    
    if saveData
    filenames{XX} = filename;
    end
    pause(1)
end

if singles || doubles
    timee = clock;
    mergeName = ['merge-', date, '-', num2str(timee(4)), '-',...
            num2str(timee(5)), '-', num2str(floor(timee(6))), '.mat'];
    mergeFiles(filenames, mergeName, sum(incl), multiplicity, deleteFiles);
end

%prepare the statistics UI table headings
stats = output.stats;
mass = output.mass;
charge = output.charge;
colmnNames = {'total number of hits'; 'total number of shots'; 'number of coincidences'; 'number unique coin'};
for i = 1:length(mass)
    var = {['mass ' num2str(mass(i)) ' charge ' num2str(charge(i))]};
    colmnNames = [colmnNames; var];
end

%output the statistics to a UI table
set(handles.stats, 'columnName', colmnNames)
set(handles.stats, 'columnWidth', 'auto')
set(handles.stats, 'data', (stats)')

if ~output.uniqueShots
    set(handles.err, 'string', 'not all of the shot numbers are unique per each coincidence');
    set(handles.cutAndPlot, 'string', 'cut and plot');
    warning('not all of the shot numbers are unique per each coincidence')
end

%reset the button string
set(handles.cutAndPlot, 'string', 'cut and plot');

%save any values saved to handles
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes when entered data in editable cell(s) in cutAndPlotParams.
function calibParams_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to cutAndPlotParams (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%make sure 'NaN' is not displayed
params = get(handles.cutAndPlotParams, 'data');
params(isnan(params)) = 0;
set(handles.cutAndPlotParams, 'data', params);

% --- Executes when entered data in editable cell(s) in cutAndPlotParams.
function cutAndPlotParams_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to cutAndPlotParams (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%make sure 'NaN' is not displayed
params = get(handles.cutAndPlotParams, 'data');
params(isnan(params)) = 0;
set(handles.cutAndPlotParams, 'data', params);

% --- Executes when entered data in editable cell(s) in prepareParams.
function prepareParams_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to prepareParams (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%make sure 'NaN' is not displayed
params = get(handles.prepareParams, 'data');
params(isnan(params)) = 0;
set(handles.prepareParams, 'data', params);

% --- Executes when entered data in editable cell(s) in commonParams.
function commonParams_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to commonParams (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%make sure 'NaN' is not displayed
params = get(handles.commonParams, 'data');
params(isnan(params)) = 0;
set(handles.commonParams, 'data', params);

% --- Executes when entered data in editable cell(s) in commonParams.
function calibPoints_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to commonParams (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


function tofNumBins_Callback(hObject, eventdata, handles)
% hObject    handle to tofNumBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of tofNumBins as text
%        str2double(get(hObject,'String')) returns contents of tofNumBins as a double

% --- Executes during object creation, after setting all properties.
function tofNumBins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tofNumBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in numHitsEqNumPart.
function numHitsEqNumPart_Callback(hObject, eventdata, handles)
% hObject    handle to numHitsEqNumPart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of numHitsEqNumPart

% --- Executes on button press in coneV.
function coneV_Callback(hObject, eventdata, handles)
% hObject    handle to coneV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of coneV

% --- Executes on button press in plTofPlot.
function plTofPlot_Callback(hObject, eventdata, handles)
% hObject    handle to plTofPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of plTofPlot

% --- Executes on button press in plTofHist.
function plTofHist_Callback(hObject, eventdata, handles)
% hObject    handle to plTofHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of plTofHist

% --- Executes on button press in plAngleHist.
function plAngleHist_Callback(hObject, eventdata, handles)
% hObject    handle to plAngleHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of plAngleHist

% --- Executes on button press in plEnergyPartHist.
function plEnergyPartHist_Callback(hObject, eventdata, handles)
% hObject    handle to plEnergyPartHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of plEnergyPartHist

% --- Executes on button press in plEnergyHist.
function plEnergyHist_Callback(hObject, eventdata, handles)
% hObject    handle to plEnergyHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of plEnergyHist

% --- Executes on button press in conHeavyPartAngle.
function conHeavyPartAngle_Callback(hObject, eventdata, handles)
% hObject    handle to conHeavyPartAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of conHeavyPartAngle

function maxHeavyPartAngle_Callback(hObject, eventdata, handles)
% hObject    handle to maxHeavyPartAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of maxHeavyPartAngle as text
%        str2double(get(hObject,'String')) returns contents of maxHeavyPartAngle as a double

% --- Executes during object creation, after setting all properties.
function maxHeavyPartAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxHeavyPartAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numCutPlotBins_Callback(hObject, eventdata, handles)
% hObject    handle to numCutPlotBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of numCutPlotBins as text
%        str2double(get(hObject,'String')) returns contents of numCutPlotBins as a double

% --- Executes during object creation, after setting all properties.
function numCutPlotBins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numCutPlotBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function tweekParams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tweekParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double

% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numExtraHits_Callback(hObject, eventdata, handles)
% hObject    handle to numExtraHits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of numExtraHits as text
%        str2double(get(hObject,'String')) returns contents of numExtraHits as a double

% --- Executes during object creation, after setting all properties.
function numExtraHits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numExtraHits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in saveData.
function saveData_Callback(hObject, eventdata, handles)
% hObject    handle to saveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of saveData

% --- Executes on button press in saveCalib.
function saveCalib_Callback(hObject, eventdata, handles)
% hObject    handle to saveCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of saveCalib

% --- Executes when selected cell(s) is changed in cutAndPlotParams.
function cutAndPlotParams_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to cutAndPlotParams (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in conDeltat.
function conDeltat_Callback(hObject, eventdata, handles)
% hObject    handle to conDeltat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of conDeltat

function deltat_Callback(hObject, eventdata, handles)
% hObject    handle to deltat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of deltat as text
%        str2double(get(hObject,'String')) returns contents of deltat as a double

% --- Executes during object creation, after setting all properties.
function deltat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function deltax_Callback(hObject, eventdata, handles)
% hObject    handle to deltax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of deltax as text
%        str2double(get(hObject,'String')) returns contents of deltax as a double

% --- Executes during object creation, after setting all properties.
function deltax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function deltay_Callback(hObject, eventdata, handles)
% hObject    handle to deltay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of deltay as text
%        str2double(get(hObject,'String')) returns contents of deltay as a double

% --- Executes during object creation, after setting all properties.
function deltay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plMomDir.
function plMomDir_Callback(hObject, eventdata, handles)
% hObject    handle to plMomDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of plMomDir

% --- Executes on button press in singles.
function singles_Callback(hObject, eventdata, handles)
% hObject    handle to singles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of singles

% --- Executes on button press in doubles.
function doubles_Callback(hObject, eventdata, handles)
% hObject    handle to doubles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of doubles

% --- Executes on button press in deleteFiles.
function deleteFiles_Callback(hObject, eventdata, handles)
% hObject    handle to deleteFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of deleteFiles

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over loadfile.
function loadfile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over findTof.
function findTof_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to findTof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on findTof and none of its controls.
function findTof_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to findTof (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on loadfile and none of its controls.
function loadfile_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on filename and none of its controls.
function filename_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function filename_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over filename.
function filename_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in shutter1.
function shutter1_Callback(hObject, eventdata, handles)
% hObject    handle to shutter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of shutter1

% --- Executes on selection change in paramChoice.
function paramChoice_Callback(hObject, eventdata, handles)
% hObject    handle to paramChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns paramChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from paramChoice

% --- Executes on selection change in shutterChoice.
function shutterChoice_Callback(hObject, eventdata, handles)
% hObject    handle to shutterChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns shutterChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from shutterChoice

% --- Executes during object creation, after setting all properties.
function shutterChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shutterChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savePrepare.
function savePrepare_Callback(hObject, eventdata, handles)
% hObject    handle to savePrepare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savePrepare


% --- Executes on button press in loadPrepare.
function loadPrepare_Callback(hObject, eventdata, handles)
% hObject    handle to loadPrepare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function filenamePrepared_Callback(hObject, eventdata, handles)
% hObject    handle to filenamePrepared (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filenamePrepared as text
%        str2double(get(hObject,'String')) returns contents of filenamePrepared as a double


% --- Executes during object creation, after setting all properties.
function filenamePrepared_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filenamePrepared (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in includePrepared.
function includePrepared_Callback(hObject, eventdata, handles)
% hObject    handle to includePrepared (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of includePrepared


% --- Executes on button press in prepTimeEst.
function prepTimeEst_Callback(hObject, eventdata, handles)
% hObject    handle to prepTimeEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of prepTimeEst


% --- Executes on button press in loadSave.
function loadSave_Callback(hObject, eventdata, handles)
% hObject    handle to loadSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadSave


% --- Executes on button press in loadCalib.
function loadCalib_Callback(hObject, eventdata, handles)
% hObject    handle to loadCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadCalib


% --- Executes on button press in breakUpData.
function breakUpData_Callback(hObject, eventdata, handles)
% hObject    handle to breakUpData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of breakUpData



function numEventTags_Callback(hObject, eventdata, handles)
% hObject    handle to numEventTags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numEventTags as text
%        str2double(get(hObject,'String')) returns contents of numEventTags as a double


% --- Executes during object creation, after setting all properties.
function numEventTags_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numEventTags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useBrokeData.
function useBrokeData_Callback(hObject, eventdata, handles)
% hObject    handle to useBrokeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useBrokeData


% --- Executes during object creation, after setting all properties.
function intensityChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intensityChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in delayChoice.
function delayChoice_Callback(hObject, eventdata, handles)
% hObject    handle to delayChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns delayChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from delayChoice


% --- Executes during object creation, after setting all properties.
function delayChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delayChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in polarChoice.
function polarChoice_Callback(hObject, eventdata, handles)
% hObject    handle to polarChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns polarChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from polarChoice


% --- Executes during object creation, after setting all properties.
function polarChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to polarChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in intensityChoice.
function intensityChoice_Callback(hObject, eventdata, handles)
% hObject    handle to intensityChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns intensityChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from intensityChoice



function chunksStartNum_Callback(hObject, eventdata, handles)
% hObject    handle to chunksStartNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chunksStartNum as text
%        str2double(get(hObject,'String')) returns contents of chunksStartNum as a double


% --- Executes during object creation, after setting all properties.
function chunksStartNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chunksStartNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autoCalib.
function autoCalib_Callback(hObject, eventdata, handles)
% hObject    handle to autoCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoCalib
