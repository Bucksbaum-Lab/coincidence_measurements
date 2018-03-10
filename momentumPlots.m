%clear all

f = momentumPlotsFunctions();

filename= ['G:\\2018_02_09\\analysis\\Acetylene_266_800_1300_4p3e9torr_80fs_DAn_DAn-cut-23-Feb-2018-22-4-41_M=1,12,13_Z=1(10eV, dpx=3)_' ...
           'OPEN.mat'];
load(filename);

[parallel_proj, perpendicular_proj]  =... 
    f.momProject(output.momXOut, output.momYOut, output.momZOut, output.mass, 1, [2, 3]);


%%
figure
f.momHists(output)

%%
figure
f.particleKERhist(1, output.partEnergyOut, output.mass, 35, 10)

%%
figure
[Xp, Yp, N_OPEN] = ...
    f.momentum2dDistPolar(parallel_proj, perpendicular_proj, [200, 0, 5], [12, 0, pi/2],...
                              'control on', 'H^{+}', 'C^{+}', 'CH^{+}');
%%
figure
f.momentum2dDistCartesian(parallel_proj, perpendicular_proj, 14, ...
                          'control off', 'H^{+}', 'C^{+}', 'CH^{+}')
%%
figure
plot((pi/48:pi/24:(pi/2-pi/48))*180/pi, sum(N_CLOSED(51:100, 1:12), 1) )
hold on
plot((pi/48:pi/24:(pi/2-pi/48))*180/pi, sum(N_OPEN((51:100), 1:12), 1))
title('Integrated from |p| = 2.5 to |p| = 3.5')
xlabel('Angle (degrees)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

%%
pl = pcolor(Xp, Yp, (N_OPEN - N_CLOSED)));
set(pl, 'EdgeColor', 'none');
pl;
axis equal tight
colormap bluewhitered(100)
title('on - off, relative')
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')

%%
