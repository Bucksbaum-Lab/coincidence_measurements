%{
%--------------------------------prepare----------------------------------%
%-------------------------------Version 4.0-------------------------------%

prepare generates matrices that will be used by cutandplot.  It generates
these matricies from position and tof data from LCLS data.  It will
generate matrices for X, Y, and Z with column representing the momentum
assuming a certain input mass.  It will also generate a vector that
indicates which hit and shot each entry comes from.

All of this is nicely opperated via the 'main.m' file which is a GUI.

%---------------------Last Edited on November 1, 2017---------------------%
%------Last Edited by and Original Author - Chelsea Liekhus-Schmaltz------% 
%}

function [EV, mom_tof, mom_x, mom_y, hitNo, shotNo, numHits, tof, rX, rY, shutterStatus,...
    intensityStatus, polarizationStatus, paramStatus, delayStatus] =...
    prepare(V1, VM, ss, t0, x0, y0, mass, charge, maxEV, tof, rX, rY, numHits, eVArray,...
    thetaArray, tof_Sim, r_Sim, shutterStatus, intensityStatus, polarizationStatus,...
    paramStatus, delayStatus, shotZero)

tof = tof-t0;
rX = rX-x0;
rY = rY-y0;

[numshots, maxions] = size(tof);

for nn = 1:length(mass)

    
    evalc('Sim = Flym_Sim(charge(nn), mass(nn), maxEV(nn), 0, 0, ss, V1, VM);');
    tof_Sim_min(nn) = reshape(Sim(2:2:end, 2), [1, 1])*10^3;
    
    evalc('Sim = Flym_Sim(charge(nn), mass(nn), maxEV(nn), 180, 0, ss, V1, VM);');
    tof_Sim_max(nn) = reshape(Sim(2:2:end, 2), [1, 1])*10^3;
    
end

tof_sensible = false([size(tof,1),1]);

for nn = 1:size(tof,2)
    
    for mm = 1:length(mass)
        
        tof_sensible = tof_sensible | ((tof_Sim_min(mm) <= tof(:, nn)) & (tof(:, nn) <= tof_Sim_max(mm)));
    end

end
    
tof = tof(tof_sensible, :);
rX = rX(tof_sensible, :);
rY = rY(tof_sensible, :);
numHits = numHits(tof_sensible, :);
shutterStatus = shutterStatus(tof_sensible, :);
intensityStatus = intensityStatus(tof_sensible, :);
polarizationStatus = polarizationStatus(tof_sensible, :);
paramStatus = paramStatus(tof_sensible, :);
delayStatus = delayStatus(tof_sensible, :);

%create the vectors that will hold the hit number and shot number
%associated with each event
hitNo = repmat(linspace(1, maxions, maxions), numshots, 1);
hitNo = (hitNo(tof_sensible, :))';
hitNo = hitNo(:);
shotNo = repmat((linspace(1, numshots, numshots))', 1, maxions)+shotZero;
shotNo = (shotNo(tof_sensible, :))';
shotNo = shotNo(:);

%create the vector that will record how many hits occur each shot

numHits = repmat(numHits, 1, maxions);
numHits = (numHits)';
numHits = numHits(:);

shutterStatus = repmat(shutterStatus, 1, maxions);
shutterStatus = (shutterStatus)';
shutterStatus = shutterStatus(:);

intensityStatus = repmat(intensityStatus, 1, maxions);
intensityStatus = (intensityStatus)';
intensityStatus = intensityStatus(:);

polarizationStatus = repmat(polarizationStatus, 1, maxions);
polarizationStatus = (polarizationStatus)';
polarizationStatus = polarizationStatus(:);

paramStatus = repmat(paramStatus, 1, maxions);
paramStatus = (paramStatus)';
paramStatus = paramStatus(:);

delayStatus = repmat(delayStatus, 1, maxions);
delayStatus = (delayStatus)';
delayStatus = delayStatus(:);

tof = tof';
rX = rX';
rY = rY';

tof = tof(:);
rX = rX(:);
rY = rY(:);

%eliminte nonsense events

cond = hitNo <= numHits;

rX = rX(cond);
rY = rY(cond);
hitNo = hitNo(cond);
shotNo = shotNo(cond);
numHits = numHits(cond);
shutterStatus = shutterStatus(cond);
intensityStatus = intensityStatus(cond);
polarizationStatus = polarizationStatus(cond);
paramStatus = paramStatus(cond);
delayStatus = delayStatus(cond);
tof = tof(cond);

%initialize the momentum matrices
EV = zeros(length(rX), length(mass));
mom_tof = zeros(length(rX), length(mass));
mom_x = zeros(length(rX), length(mass));
mom_y = zeros(length(rX), length(mass));

%find the momentum associated with each event assuming it has a given mass
%determined by which column it is in
'make sure to return this to a parfor: prepare line 128'
for ii = 1:length(mass)
[EV(:, ii), mom_tof(:, ii), mom_x(:, ii), mom_y(:, ii)] =...
        convertToEnergy(tof, rX, rY, eVArray(:, ii), thetaArray(:, ii), tof_Sim(:, :, ii), r_Sim(:, :, ii), mass(ii));
end
