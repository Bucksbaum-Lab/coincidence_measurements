function [eVx, eVy, eVtof] = ytof2eVy(tof, rX, rY, V1, VM, ss, charge, mass)

eVxArray = [logspace(-20, -26, 1000), -1*logspace(-26, -19, 1000)];
eVyArray = [logspace(-20, -26, 1000), -1*logspace(-26, -19, 1000)];
eVtofArray = [logspace(-20, -26, 1000), -1*logspace(-26, -19, 1000)];
[rXArray, rYArray, tofArray] = calculatedtime(V1, VM, ss, A, L, C, mass, charge, eVxArray, eVyArray, eVtofArray);

%figure()
%plot(tofArray, momZArray, '.', 'markersize', 1)
eVy = spline(tofArray, eVtofArray, tof);