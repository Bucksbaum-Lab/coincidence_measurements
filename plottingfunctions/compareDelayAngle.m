function compareDelayAngle(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, intensity)

for nn = 1:length(plotting_data)
    for mm = nn+1:length(plotting_data)
        
        if plotting_data(nn).shutter == plotting_data(mm).shutter
            
            polarDist = plotting_data(nn).polarDistNorm;
            polarDist(:,end) = [];
            
            figure
            %errorbar(linspace(0,pi,numAngle), sum(polarDist, 1), sqrt(sum(polarDist, 1)) )
            plot(linspace(0,pi,numAngle), sum(polarDist, 1))
            hold on
            
            polarDist = plotting_data(mm).polarDistNorm;
            polarDist(:,end) = [];
            
            %errorbar(linspace(0,pi,numAngle), sum(polarDist, 1), sqrt(sum(polarDist, 1)) )
            plot(linspace(0,pi,numAngle), sum(polarDist, 1))
            
            title(['Integrated from |p| = 0 to |p| = ' num2str(plotting_data(nn).rs(end)) ', '...
                char(shutterTitle(strcmp(plotting_data(nn).shutter,shutterMarkers))) ', intensity ' intensity])
            xlabel('Angle (rad)')
            ylabel('counts')
            legend(char(delayTitle(strcmp(plotting_data(nn).delay,delayMarkers))),...
                char(delayTitle(strcmp(plotting_data(mm).delay,delayMarkers))),'Location','southeast')
            grid on
        
        end
    end
end
