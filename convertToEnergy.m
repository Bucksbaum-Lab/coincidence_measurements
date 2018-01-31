function [EV, mom_tof, mom_x, mom_y, Theta] = convertToEnergy(tof, rX, rY, V1, VM, ss, charge, mass, maxEV, EVlength, Thetalength)

eVArray = linspace(0, maxEV, EVlength);

thetaArray = linspace(0, 180, Thetalength);

%{
%estimate how long it will take to load file
estimatedTime = 0;

for ii = 1:10
    tic
    evalc('Sim = Flym_Sim(1, [12,12,12,12,12,12,12,12,12,12,12], 0, 0, 0, ss, V1, VM);');
    estimatedTime = estimatedTime+toc;
    tic
    evalc('Sim = Flym_Sim(1, 12, 0, 0, 0, ss, V1, VM);');
    estimatedTime = estimatedTime-toc;
end

estimatedTime = estimatedTime/100*length(eVArray)*length(thetaArray)*2.5;
fclose('all');
clear testLoad;

['estimated completed simulate time: ' datestr(datetime('now')+seconds(estimatedTime))]
%}

%do the simulations
evalc('Sim = Flym_Sim(charge, mass, eVArray, thetaArray, 0, ss, V1, VM);');

tof_Sim = reshape(Sim(2:2:end, 2), [length(eVArray), length(thetaArray)])*10^3;
r_Sim = abs(reshape(Sim(2:2:end, 3)-55, [length(eVArray), length(thetaArray)]));

%find the EV and theta
rR = sqrt(rX.^2+rY.^2);

[thetaArray, eVArray] = meshgrid(thetaArray, eVArray);

EV = griddata(tof_Sim, r_Sim, eVArray, tof, rR, 'v4');
Theta = griddata(tof_Sim, r_Sim, thetaArray, tof, rR, 'v4')*pi()/180;


%get rid of nonsense answers
nonSense_tof_max = max(max(tof_Sim));
nonSense_tof_min = min(min(tof_Sim));
EV((tof > nonSense_tof_max)|(tof < nonSense_tof_min)) = NaN;
Theta((tof > nonSense_tof_max)|(tof < nonSense_tof_min)) = NaN;
%}

%find the momentum vectors as well
MomTotal = sqrt(2*mass*EV);
mom_tof = MomTotal.*cos(Theta);
mom_x = MomTotal.*sin(Theta).*rX./rR;
mom_y = MomTotal.*sin(Theta).*rY./rR;
