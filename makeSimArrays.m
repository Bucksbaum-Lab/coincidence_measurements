function [eVArray, thetaArray, tof_Sim, r_Sim] =...
    makeSimArrays(V1, VM, ss, mass, charge, maxEV, EVlength, Thetalength)

for mm=1:length(mass)

    eVArray(:,mm) = (linspace(0, (maxEV(mm))^.5, EVlength)).^2;
    %eVArray(:,mm) = (linspace(0, (maxEV(mm)).^2, EVlength)).^.5;
    %eVArray(:,mm) = linspace(0, maxEV(mm), EVlength);

    thetaArray(:,mm) = linspace(0, 180, Thetalength);
    %thetaArray(:,mm) = (linspace(0, sqrt(180), Thetalength)).^2;
    %thetaArray(:,mm) = linspace(0, 180, Thetalength);

    %do the simulations
    evalc('Sim = Flym_Sim(charge(mm), mass(mm), eVArray(:,mm), thetaArray(:,mm), 0, ss, V1, VM);');

    tof_Sim(:,:,mm) = reshape(Sim(2:2:end, 2), [EVlength, Thetalength])*10^3;
    r_Sim(:,:,mm) = abs(reshape(Sim(2:2:end, 4)-55, [EVlength, Thetalength]));
    
end