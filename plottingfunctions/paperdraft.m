%clear all

f = momentumPlotsFunctions();

delay = 1000;
intensity = 'high';

filename_closed_pd= ['G:\2018_05_03\analysis\'...
    'acetelene_1ps_6p3n10torr_0p6nd_DAn_DAn-loaded_data-10-May-2018-22-10-25_-cut-closed_' intensity '_1000fs.mat'];
load(filename_closed_pd);

output_closed_pd = output;

filename_open_pd= ['G:\2018_05_03\analysis\'...
    'acetelene_1ps_6p3n10torr_0p6nd_DAn_DAn-loaded_data-10-May-2018-22-10-25_-cut-open_' intensity '_1000fs.mat'];
load(filename_open_pd);

output_open_pd = output;

filename_closed_nd= ['G:\2018_05_03\analysis\'...
    'acetelene_1ps_6p3n10torr_0p6nd_DAn_DAn-loaded_data-10-May-2018-22-10-25_-cut-closed_' intensity '_n1000fs.mat'];
load(filename_closed_nd);

output_closed_nd = output;

filename_open_nd= ['G:\2018_05_03\analysis\'...
    'acetelene_1ps_6p3n10torr_0p6nd_DAn_DAn-loaded_data-10-May-2018-22-10-25_-cut-open_' intensity '_n1000fs.mat'];
load(filename_open_nd);

output_open_nd = output;

clear('output')

[parallel_proj_closed_pd, perpendicular_proj_closed_pd]  =... 
    f.momProject(output_closed_pd.momXOut, output_closed_pd.momYOut, output_closed_pd.momZOut, output_closed_pd.mass, 1, [2, 3]);

[parallel_proj_open_pd, perpendicular_proj_open_pd]  =... 
    f.momProject(output_open_pd.momXOut, output_open_pd.momYOut, output_open_pd.momZOut, output_open_pd.mass, 1, [2, 3]);

[parallel_proj_closed_nd, perpendicular_proj_closed_nd]  =... 
    f.momProject(output_closed_nd.momXOut, output_closed_nd.momYOut, output_closed_nd.momZOut, output_closed_nd.mass, 1, [2, 3]);

[parallel_proj_open_nd, perpendicular_proj_open_nd]  =... 
    f.momProject(output_open_nd.momXOut, output_open_nd.momYOut, output_open_nd.momZOut, output_open_nd.mass, 1, [2, 3]);

%%
figure
f.momHists(output_closed_pd)
mtit(['momentum sums control off ' num2str(delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_closed_pd))])

figure
f.momHists(output_open_pd)
mtit(['momentum sums control on ' num2str(delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_open_pd))])

figure
f.momHists(output_closed_nd)
mtit(['momentum sums control off ' num2str(-delay) ' fs; intensity ' intensity ' total counts ' num2str(length(parallel_proj_closed_nd))])

figure
f.momHists(output_open_nd)
mtit(['momentum sums control on ' num2str(-delay) ' fs; intensity ' intensity ' total counts ' num2str(length(parallel_proj_open_nd))])

%%
figure
f.particleKERhist(1, output_closed_pd.partEnergyOut, output_closed_pd.mass, 30, 30,...
    ['control off ' num2str(delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_closed_pd))])

figure
f.particleKERhist(1, output_open_pd.partEnergyOut, output_open_pd.mass, 30, 30,...
    ['control on ' num2str(delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_open_pd))])

figure
f.particleKERhist(1, output_closed_nd.partEnergyOut, output_closed_nd.mass, 30, 30,...
    ['control off ' num2str(-delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_closed_nd))])

figure
f.particleKERhist(1, output_open_nd.partEnergyOut, output_open_nd.mass, 30, 30,...
    ['control on ' num2str(-delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_open_nd))])

%%

numMom = 5;
maxMom = 8;
numAngle = 5;

figure
[Xpd, Ypd, N_CLOSED_pd, rs] = ...
    f.momentum2dDistPolar(parallel_proj_closed_pd, perpendicular_proj_closed_pd, [numMom, 0, maxMom], [numAngle, 0, pi/2],...
    ['control off ' num2str(delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_closed_pd))],...
    'H^{+}', 'C^{+}', 'CH^{+}');
colorbar                          

figure
[~, ~, N_OPEN_pd, ~] = ...
    f.momentum2dDistPolar(parallel_proj_open_pd, perpendicular_proj_open_pd, [numMom, 0, maxMom], [numAngle, 0, pi/2],...
    ['control on ' num2str(delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_open_pd))],...
'H^{+}', 'C^{+}', 'CH^{+}');
colorbar

figure
[Xnd, Ynd, N_CLOSED_nd, rs] = ...
    f.momentum2dDistPolar(parallel_proj_closed_nd, perpendicular_proj_closed_nd, [numMom, 0, maxMom], [numAngle, 0, pi/2],...
    ['control off ' num2str(-delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_closed_nd))],...
    'H^{+}', 'C^{+}', 'CH^{+}');
colorbar                          

figure
[~, ~, N_OPEN_nd, ~] = ...
    f.momentum2dDistPolar(parallel_proj_open_nd, perpendicular_proj_open_nd, [numMom, 0, maxMom], [numAngle, 0, pi/2],...
    ['control on ' num2str(-delay) ' fs; intensity ' intensity '; total counts ' num2str(length(parallel_proj_open_nd))],...
    'H^{+}', 'C^{+}', 'CH^{+}');
colorbar

%%
figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(round(3/maxMom*numMom))) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = ' num2str(rs(round(3/maxMom*numMom))) ' to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_pd, 1), sqrt(sum(N_CLOSED_pd, 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_pd, 1), sqrt(sum(N_OPEN_pd, 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(round(3/maxMom*numMom))) ', ' num2str(-delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = ' num2str(rs(round(3/maxMom*numMom))) ' to |p| = ' num2str(rs(end)) ', ' num2str(-delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_nd, 1), sqrt(sum(N_CLOSED_nd, 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_nd, 1), sqrt(sum(N_OPEN_nd, 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(end)) ', ' num2str(-delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

%%
figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(round(3/maxMom*numMom))) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control on, positive delay','control on, negative delay','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_OPEN_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = ' num2str(rs(round(3/maxMom*numMom))) ' to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control on, positive delay','control on, negative delay','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_pd, 1), sqrt(sum(N_OPEN_pd, 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN_nd, 1), sqrt(sum(N_OPEN_nd, 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control on, positive delay','control on, negative delay','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_pd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_nd(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(round(3/maxMom*numMom))) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off, positive delay','control off, negative delay','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_pd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED_nd(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = ' num2str(rs(round(3/maxMom*numMom))) ' to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off, positive delay','control off, negative delay','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_pd, 1), sqrt(sum(N_CLOSED_pd, 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED_nd, 1), sqrt(sum(N_CLOSED_nd, 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs; intensity ' intensity])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off, positive delay','control off, negative delay','Location','southeast')
grid on

%%
figure
pl = pcolor(Xpd, Ypd, (N_OPEN_pd/sum(sum(N_OPEN_pd)) - N_CLOSED_pd/sum(sum(N_CLOSED_pd))));
set(pl, 'EdgeColor', 'none');
pl;
axis equal tight
colormap bluewhitered(100)
title(['on - off, relative, ' num2str(delay) ' fs; intensity ' intensity])
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')

figure
pl = pcolor(Xpd, Ypd, (N_OPEN_nd/sum(sum(N_OPEN_nd)) - N_CLOSED_pd/sum(sum(N_CLOSED_nd))));
set(pl, 'EdgeColor', 'none');
pl;
axis equal tight
colormap bluewhitered(100)
title(['on - off, relative, ' num2str(-delay) ' fs; intensity ' intensity])
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')

%%
caxislimit = .025;
figure
pl = pcolor(Xpd, Ypd, (N_OPEN_pd/sum(sum(N_OPEN_pd)) - N_OPEN_nd/sum(sum(N_OPEN_nd))));
set(pl, 'EdgeColor', 'none');
pl;
axis equal tight
colormap bluewhitered(100)
title(['positive - negative delay, control on; intensity ' intensity])
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
colorbar
caxis([-caxislimit, caxislimit])

figure
pl = pcolor(Xpd, Ypd, (N_CLOSED_pd/sum(sum(N_CLOSED_pd)) - N_CLOSED_nd/sum(sum(N_CLOSED_nd))));
set(pl, 'EdgeColor', 'none');
pl;
axis equal tight
colormap bluewhitered(100)
title(['positive - negative delay, control off; intensity ' intensity])
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
colorbar
caxis([-caxislimit, caxislimit])
%%
