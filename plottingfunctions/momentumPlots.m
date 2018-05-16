%clear all

f = momentumPlotsFunctions();

filename= ['C:\Data\2018_03_09\analysis\'...
    'Acetylene_5p4en10torr_266_5p6mW_1300_38mW_800_103mW_500fsDelay_DAn_DAn-loaded_data-18-Apr-2018-16-18-29_-cut-open.mat'];
load(filename);

[parallel_proj, perpendicular_proj]  =... 
    f.momProject(output.momXOut, output.momYOut, output.momZOut, output.mass, 1, [2, 3]);


%%
figure
f.momHists(output)

%%
figure
f.particleKERhist(1, output.partEnergyOut, output.mass, 60, 30)

%%
figure
[Xp, Yp, N_OPEN] = ...
    f.momentum2dDistPolar(parallel_proj, perpendicular_proj, [15, 0, 8], [6, 0, pi/2],...
                              'control on', 'H^{+}', 'C^{+}', 'CH^{+}');
%%
figure
f.momentum2dDistCartesian(parallel_proj, perpendicular_proj, 14, ...
                          'control off', 'H^{+}', 'C^{+}', 'CH^{+}')
%%
figure
%(pi/48:pi/24:(pi/2-pi/48))*180/pi, 
errorbar(sum(N_CLOSED(5:5, 1:7), 1), 1./sqrt(sum(N_CLOSED(5:5, 1:7), 1)) )
hold on
errorbar(sum(N_OPEN(5:5, 1:7), 1), 1./sqrt(sum(N_OPEN(5:5, 1:7), 1)) )
title('Integrated from |p| = 2.13 to |p| = 4.26')
xlabel('Angle (degrees)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

%%
figure
pl = pcolor(Xp, Yp, (N_OPEN/sum(sum(N_OPEN)) - N_CLOSED/sum(sum(N_CLOSED))));
set(pl, 'EdgeColor', 'none');
pl;
axis equal tight
colormap bluewhitered(100)
title('on - off, relative')
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')

%%
