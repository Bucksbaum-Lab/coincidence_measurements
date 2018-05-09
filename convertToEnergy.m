function [EV, mom_tof, mom_x, mom_y, Theta] = convertToEnergy(tof, rX, rY, eVArray, thetaArray, tof_Sim, r_Sim, mass)

%find the EV and theta
rR = sqrt(rX.^2+rY.^2);

[thetaArray, eVArray] = meshgrid(thetaArray, eVArray);

nonSense_tof_max = max(max(tof_Sim));
nonSense_tof_min = min(min(tof_Sim));
isresonable_tof = (nonSense_tof_min <= tof)&(tof <= nonSense_tof_max);

EV = nan(size(tof));
Theta = nan(size(tof));

EV(isresonable_tof)    = griddata(tof_Sim, r_Sim, eVArray,...
                                  tof(isresonable_tof), rR(isresonable_tof), 'linear');
Theta(isresonable_tof) = griddata(tof_Sim, r_Sim, thetaArray,...
                                  tof(isresonable_tof), rR(isresonable_tof), 'linear')*pi()/180;

%find the momentum vectors as well
MomTotal = sqrt(2*mass*EV);
mom_tof = MomTotal.*cos(Theta);
mom_x = MomTotal.*sin(Theta).*rX./rR;
mom_y = MomTotal.*sin(Theta).*rY./rR;

