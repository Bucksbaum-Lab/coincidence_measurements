function plotting_data = anglesAndProjections(plotting_data)

for nn = 1:length(plotting_data)

    [plotting_data(nn).parallel_proj, plotting_data(nn).perpendicular_proj]  =... 
        momProject(plotting_data(nn).output.momXOut, plotting_data(nn).output.momYOut, plotting_data(nn).output.momZOut, plotting_data(nn).output.mass, 1, [2, 3]);

    plotting_data(nn).theta = atan2(plotting_data(nn).perpendicular_proj, plotting_data(nn).parallel_proj);
    
end
