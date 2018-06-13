function compareShutterAngle(plotting_data, delayMarkers, shutterMarkers, delayTitle, shutterTitle, numAngle, intensity)

for nn = 1:length(plotting_data)
    for mm = nn+1:length(plotting_data)
        
        if plotting_data(nn).delay == plotting_data(mm).delay
            
            polarDist = plotting_data(nn).polarDist;
            polarDist(:,end) = [];
            
            figure
            errorbar(linspace(0,pi,numAngle), sum(polarDist, 1), sqrt(sum(polarDist, 1)) )
            hold on
            
            polarDist = plotting_data(mm).polarDist;
            polarDist(:,end) = [];
            
            errorbar(linspace(0,pi,numAngle), sum(polarDist, 1), sqrt(sum(polarDist, 1)) )
            
            title(['Integrated from |p| = 0 to |p| = ' num2str(plotting_data(nn).rs(end)) ', '...
                char(delayTitle(strcmp(plotting_data(nn).delay,delayMarkers))) ', intensity ' intensity])
            xlabel('Angle (rad)')
            ylabel('counts')
            legend(char(shutterTitle(strcmp(plotting_data(nn).shutter,shutterMarkers))),...
                char(shutterTitle(strcmp(plotting_data(mm).shutter,shutterMarkers))),'Location','southeast')
            grid on
        
        end
    end
end
