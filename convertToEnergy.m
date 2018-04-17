function [EV, mom_tof, mom_x, mom_y, Theta] = ...
    convertToEnergy(tof, rX, rY, V1, VM, ss, charge, mass, maxEV, EVlength, Thetalength, uniqueid)

eVArray = linspace(0, maxEV, EVlength);

thetaArray = linspace(0, 180, Thetalength);

%do the simulations
evalc('Sim = Flym_Sim(charge, mass, eVArray, thetaArray, 0, ss, V1, VM, uniqueid);');

tof_Sim = reshape(Sim(2:2:end, 2), [length(eVArray), length(thetaArray)])*10^3;
r_Sim = abs(reshape(Sim(2:2:end, 3)-55, [length(eVArray), length(thetaArray)]));

%find the EV and theta
rR = sqrt(rX.^2+rY.^2);

[thetaArray, eVArray] = meshgrid(thetaArray, eVArray);

nonSense_tof_max = max(max(tof_Sim));
nonSense_tof_min = min(min(tof_Sim));
isresonable_tof = (nonSense_tof_min <= tof)&(tof <= nonSense_tof_max);

EV = nan(size(tof));
Theta = nan(size(tof));

EV(isresonable_tof)    = griddata(tof_Sim, r_Sim, eVArray,...
                                  tof(isresonable_tof), rR(isresonable_tof), 'v4');
Theta(isresonable_tof) = griddata(tof_Sim, r_Sim, thetaArray,...
                                  tof(isresonable_tof), rR(isresonable_tof), 'v4')*pi()/180;

%find the momentum vectors as well
MomTotal = sqrt(2*mass*EV);
mom_tof = MomTotal.*cos(Theta);
mom_x = MomTotal.*sin(Theta).*rX./rR;
mom_y = MomTotal.*sin(Theta).*rY./rR;

