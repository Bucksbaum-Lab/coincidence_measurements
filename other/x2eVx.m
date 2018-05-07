function momZ = tof2momz(tof, V1, VM, s, A, L, C, charge, mass)

momZArray = [logspace(-20, -26, 1000), -1*logspace(-26, -19, 1000)];
tofArray = calculatedtime(V1, VM, s, A, L, C, mass, momZArray, charge);

%figure()
%plot(tofArray, momZArray, '.', 'markersize', 1)
momZ = spline(tofArray, momZArray, tof);