function plotting_data = polarPlots(plotting_data, numMom, maxMom, numAngle, shutterTitle, delayTitle, intensity)

ll = 0;

for nn = 1:length(shutterTitle)
    for mm = 1:length(delayTitle)

        ll = ll+1;
        
        figure
        [plotting_data(ll).Xpolar, plotting_data(ll).Ypolar, plotting_data(ll).polarDist, plotting_data(ll).rs] = ...
            momentum2dDistPolar(plotting_data(ll).parallel_proj, plotting_data(ll).perpendicular_proj,...
            [numMom, 0, maxMom], [numAngle, 0, pi],...
            [char(shutterTitle(nn)) ', ' char(delayTitle(mm)) ', intensity ' intensity...
            ', total counts ' num2str(length(plotting_data(ll).parallel_proj))],...
            'H^{+}', 'C^{+}', 'CH^{+}');
        colorbar
    
    end
end
