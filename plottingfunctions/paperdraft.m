%clear all

f = momentumPlotsFunctions();

delay = 500;

filename_closed= ['C:\Data\2018_03_09\analysis\'...
    'Acetylene_5p4en10torr_266_5p6mW_1300_38mW_800_103mW_500fsDelay_DAn_DAn-loaded_data-18-Apr-2018-16-18-29_-cut-closed.mat'];
load(filename_closed);

output_closed = output;

filename_open= ['C:\Data\2018_03_09\analysis\'...
    'Acetylene_5p4en10torr_266_5p6mW_1300_38mW_800_103mW_500fsDelay_DAn_DAn-loaded_data-18-Apr-2018-16-18-29_-cut-open.mat'];
load(filename_open);

output_open = output;

clear('output')

[parallel_proj_closed, perpendicular_proj_closed]  =... 
    f.momProject(output_closed.momXOut, output_closed.momYOut, output_closed.momZOut, output_closed.mass, 1, [2, 3]);

[parallel_proj_open, perpendicular_proj_open]  =... 
    f.momProject(output_open.momXOut, output_open.momYOut, output_open.momZOut, output_open.mass, 1, [2, 3]);

%%
figure
f.momHists(output_closed)
mtit(['momentum sums control off ' num2str(delay) ' fs'])

figure
f.momHists(output_open)
mtit(['momentum sums control on ' num2str(delay) ' fs'])

%%
figure
f.particleKERhist(1, output_closed.partEnergyOut, output_closed.mass, 30, 30, ['control off ' num2str(delay) ' fs'])

figure
f.particleKERhist(1, output_open.partEnergyOut, output_open.mass, 30, 30, ['control on ' num2str(delay) ' fs'])

%%

numMom = 10;
maxMom = 8;
numAngle = 6;

figure
[Xp, Yp, N_CLOSED, rs] = ...
    f.momentum2dDistPolar(parallel_proj_closed, perpendicular_proj_closed, [numMom, 0, maxMom], [numAngle, 0, pi/2],...
    ['control off ' num2str(delay) ' fs'], 'H^{+}', 'C^{+}', 'CH^{+}');
colorbar                          

figure
[~, ~, N_OPEN, ~] = ...
    f.momentum2dDistPolar(parallel_proj_open, perpendicular_proj_open, [numMom, 0, maxMom], [numAngle, 0, pi/2],...
    ['control on ' num2str(delay) ' fs'], 'H^{+}', 'C^{+}', 'CH^{+}');
colorbar

%%
figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1), sqrt(sum(N_OPEN(1:round(3/maxMom*numMom-1), 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(round(3/maxMom*numMom))) ', ' num2str(delay) ' fs'])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_CLOSED(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1), sqrt(sum(N_OPEN(round(3/maxMom*numMom):end, 1:(numAngle+1)), 1)) )
title(['Integrated from |p| = ' num2str(rs(round(3/maxMom*numMom))) ' to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs'])
xlabel('Angle (rad)')
ylabel('counts')
legend('control off','control on','Location','southeast')
grid on

figure
errorbar(linspace(0,pi/2,numAngle+1), sum(N_CLOSED, 1), sqrt(sum(N_CLOSED, 1)) )
hold on
errorbar(linspace(0,pi/2,numAngle+1), sum(N_OPEN, 1), sqrt(sum(N_OPEN, 1)) )
title(['Integrated from |p| = 0 to |p| = ' num2str(rs(end)) ', ' num2str(delay) ' fs'])
xlabel('Angle (rad)')
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
title(['on - off, relative, ' num2str(delay) ' fs'])
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')

%%
