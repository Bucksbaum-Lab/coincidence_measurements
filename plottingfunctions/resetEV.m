function plotting_data = resetEV(plotting_data)

for nn = 1:length(plotting_data)

    plotting_data(nn).output.partEnergyOut = ...
        (plotting_data(nn).output.momXOut.^2+plotting_data(nn).output.momYOut.^2+...
        plotting_data(nn).output.momZOut.^2)./plotting_data(nn).output.mass/2;
    
    plotting_data(nn).output.KER = sum(plotting_data(nn).output.partEnergyOut,2);

end