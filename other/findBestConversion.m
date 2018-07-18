clear
clc

EVlength_det = 300;
thetalength_det = 200;

[eVArray_detailed, thetaArray_detailed, tof_Sim_detailed, r_Sim_detailed] =...
    makeSimArrays(4000, -2250, 37.5, 1, 1, 20, EVlength_det, thetalength_det);
%%
tic
[eVArray, thetaArray, tof_Sim, r_Sim] = makeSimArrays(4000, -2250, 37.5, 1, 1, 20, 25, 80);

thetaMatrix_detailed = repmat(thetaArray_detailed,1,EVlength_det);
thetaMatrix_detailed = thetaMatrix_detailed';
thetaMatrix_detailed = thetaMatrix_detailed(:);
eVMatrix_detailed = repmat(eVArray_detailed',thetalength_det,1);
eVMatrix_detailed = eVMatrix_detailed';
eVMatrix_detailed = eVMatrix_detailed(:);

[EV_det, mom_tof_det, mom_x_det, mom_y_det, Theta_det] =...
    convertToEnergy(tof_Sim_detailed(:), r_Sim_detailed(:),...
    zeros(size(r_Sim_detailed(:))), eVArray, thetaArray, tof_Sim, r_Sim, 1);

Theta_det = Theta_det(~isnan(EV_det))*180/pi;
thetaMatrix_detailed = thetaMatrix_detailed(~isnan(EV_det));
eVMatrix_detailed = eVMatrix_detailed(~isnan(EV_det));
EV_det = EV_det(~isnan(EV_det));

Test_res_eV = (eVMatrix_detailed-EV_det);
Test_Error_eV = sum((Test_res_eV).^2)/(EVlength_det*thetalength_det);
toc
figure()
plot(eVMatrix_detailed, Test_res_eV,'o')

Test_res_t = (thetaMatrix_detailed-Theta_det);
Test_Error_t = sum((Test_res_t).^2)/(EVlength_det*thetalength_det);

figure()
plot(thetaMatrix_detailed, Test_res_t,'o')
