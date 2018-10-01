function [tof_Sim, x_Sim, y_Sim, R0pos, T0pos] =...
    makeSimArrays2(V1, VM, ss, mass, charge, deltaR0, R0length)

R0pos = linspace(-deltaR0, deltaR0, R0length) + 55;
T0pos = linspace(-deltaR0, deltaR0, R0length) + ss;

%do the simulations
evalc('Sim = Flym_Sim2(charge, mass, T0pos, R0pos, R0pos, V1, VM);');

tof_Sim = reshape(Sim(2:2:end, 2), [R0length, R0length, R0length])*10^3;
x_Sim = reshape(Sim(2:2:end, 3)-55, [R0length, R0length, R0length]);
y_Sim = reshape(Sim(2:2:end, 4)-55, [R0length, R0length, R0length]);

R0pos = R0pos - 55;
T0pos = T0pos - ss;