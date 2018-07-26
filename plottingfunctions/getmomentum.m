function [plotting_data] = getmomentum(folder, filename, intensity, delay, limit, shutterMarkers, delayMarkers)

ll = 0;

for nn = 1:length(shutterMarkers)
    for mm = 1:length(delayMarkers)
        
        ll = ll+1;

        loadname = [folder filename char(shutterMarkers(nn)) intensity char(delayMarkers(mm)) delay 'fs.mat'];
        load(loadname);

        cond = (sum(output.momXOut,2)<limit) & (sum(output.momYOut,2)<limit) & (sum(output.momZOut,2)<limit);

        output.momXOut =...
            getCOMmomentum(output.mass(1), output.momXOut(cond,1), output.mass(2), output.momXOut(cond,2), output.mass(3), output.momXOut(cond,3));

        output.momYOut =...
            getCOMmomentum(output.mass(1), output.momYOut(cond,1), output.mass(2), output.momYOut(cond,2), output.mass(3), output.momYOut(cond,3));

        output.momZOut =...
            getCOMmomentum(output.mass(1), output.momZOut(cond,1), output.mass(2), output.momZOut(cond,2), output.mass(3), output.momZOut(cond,3));

        output.hitNoOut = output.hitNoOut(cond,:);
        output.shotNoOut = output.shotNoOut(cond,:);
        output.numHitsOut = output.numHitsOut(cond,:);
        output.tofOut = output.tofOut(cond,:);
        output.partEnergyOut = output.partEnergyOut(cond,:);
        output.KER = output.KER(cond,:);
        
        plotting_data(ll).output = output;
        plotting_data(ll).shutter = shutterMarkers(nn);
        plotting_data(ll).delay = delayMarkers(mm);

    end
end

