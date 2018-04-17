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

function [EV, mom_tof, mom_x, mom_y, hitNo, shotNo, numHits, tof, rX, rY, shutter, full, low] =...
    prepare(V1, VM, ss, t0, x0, y0, mass, charge, maxEV, tof, rX, rY, numHits, eVArray,...
    thetaArray, tof_Sim, r_Sim, shutter, full, low)

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
shutter = shutter(tof_sensible, :);
full = full(tof_sensible, :);
low = low(tof_sensible, :);

%create the vectors that will hold the hit number and shot number
%associated with each event
hitNo = repmat(linspace(1, maxions, maxions), numshots, 1);
hitNo = (hitNo(tof_sensible, :))';
hitNo = hitNo(:);
shotNo = repmat((linspace(1, numshots, numshots))', 1, maxions);
shotNo = (shotNo(tof_sensible, :))';
shotNo = shotNo(:);

%create the vector that will record how many hits occur each shot

numHits = repmat(numHits, 1, maxions);
numHits = (numHits)';
numHits = numHits(:);

shutter = repmat(shutter, 1, maxions);
shutter = (shutter)';
shutter = shutter(:);

full = repmat(full, 1, maxions);
full = (full)';
full = full(:);

low = repmat(low, 1, maxions);
low = (low)';
low = low(:);

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
shutter = shutter(cond);
full = full(cond);
low = low(cond);
tof = tof(cond);
%}


%initialize the momentum matrices
EV = zeros(length(rX), length(mass));
mom_tof = zeros(length(rX), length(mass));
mom_x = zeros(length(rX), length(mass));
mom_y = zeros(length(rX), length(mass));

%find the momentum associated with each event assuming it has a given mass
%determined by which column it is in

parfor ii = 1:length(mass)
[EV(:, ii), mom_tof(:, ii), mom_x(:, ii), mom_y(:, ii)] =...
        convertToEnergy(tof, rX, rY, eVArray(:, ii), thetaArray(:, ii), tof_Sim(:, :, ii), r_Sim(:, :, ii), mass(ii));end
