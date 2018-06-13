function plotting_data = polarPlotsNormalized(plotting_data, numMom, maxMom, numAngle, shutterTitle, delayTitle, intensity)

ll = 0;

cmax = 0;

for nn = 1:length(shutterTitle)
    for mm = 1:length(delayTitle)

        ll = ll+1;
        
        figure
        [plotting_data(ll).XpolarNorm, plotting_data(ll).YpolarNorm, plotting_data(ll).polarDistNorm, plotting_data(ll).rsNorm] = ...
            momentum2dDistPolarNormalized(plotting_data(ll).parallel_proj, plotting_data(ll).perpendicular_proj,...
            [numMom, 0, maxMom], [numAngle, 0, pi],...
            [char(shutterTitle(nn)) ', ' char(delayTitle(mm)) ', intensity ' intensity...
            ', total counts ' num2str(length(plotting_data(ll).parallel_proj))],...
            'H^{+}', 'C^{+}', 'CH^{+}');
        colorbar
    
        cmax = max(cmax, max(max(plotting_data(ll).polarDistNorm)));
        
    end
end

h = findobj('type', 'figure');

for nn = 1:ll
    
    figure(h(nn))
    caxis([0 cmax])
    
end