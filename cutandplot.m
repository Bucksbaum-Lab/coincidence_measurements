%{
%-------------------------------cutandplot--------------------------------%
%-------------------------------Version 4.0-------------------------------%

cutandplot will take data from the prepare function and cut the data
according to the options listed in the input of the function.  It is 
assumed that the user is working with C2H2 or C2D2.  

The possible cuts are
1) requiring that there are only x number of hits in each shot than the
number of particles in the desired coincidence.  2) requiring that the 
total momentum associated with each coincidence is bounded. 3) requiring
that the angle between the two heaviest particles in each shot is
constrained.  4) requiring that the angle between the lightest particle(s)
and the heaviest particles is bounded.  and 5) requiring that the total
energy associated with each shot is greater than some value.

After these cuts are made, cutandplot generates plots from a selection of
plots.  

All of this is nicely opperated via the 'main.m' file which is a GUI.

%}

function [output] = cutandplot(mass, charge, colms, momX, momY, momZ, partEnergy,...
    tof, hitNo, shotNo, numHits, numExtraHits, coneV, maxeV,...
    plTofPlot, plTofHist, plEnergyPartHist, plEnergyHist, plMomDir,...
    numCutPlotBins, deltat, deltax, deltay, conDeltat, incl)

%-----------------------------------cuts----------------------------------%

%use the column vector to figure out how to get properly label the incoming
%momentum matrices
mQ = mass./charge;
numPart = length(mQ);
[~,R] = sort(mQ);
colms = colms(R);
mass = mass(R);
charge = charge(R);
incl = incl(R);
maxeV = maxeV(R);

%initialize any necessary vectors and matrices
stats = zeros(numPart+4, 1);
stats(1) = length(tof);
stats(2) = length(unique(shotNo));

%record how many of each particle there is before any cuts
for i = 5:numPart+4
    stats(i) = sum(sum(~isnan(partEnergy)));
end

%first cut: require that the total number of hits equals the total number
%of particles in the desired coincidence
[momXOut, momYOut, momZOut, partEnergyOut, tofOut,...
    hitNoOut,shotNoOut, numHitsOut, cutTof] = ...
    FirstCut(numPart, numExtraHits, colms, momX, momY, momZ, ...
             partEnergy, tof, hitNo, shotNo, numHits);

[~, R] = sort(incl, 'descend');

momXOut = momXOut(:, R);
momYOut = momYOut(:, R);
momZOut = momZOut(:, R);
tofOut = tofOut(:, R);
hitNoOut = hitNoOut(:, R);
shotNoOut = shotNoOut(:, R);
numHitsOut = numHitsOut(:, R);
partEnergyOut = partEnergyOut(:, R);
maxeV = maxeV(R);
mass = mass(R);
charge = charge(R);
incl = incl(R);

clear momX
clear momY
clear momZ
clear partEnergy

%second cut: constrain the energy of each particle by the tof
if coneV && ~(isempty(maxeV))
    cond = true(size(momXOut,1),1);
    
    for i = 1:numPart
        cond = cond & (partEnergyOut(:, i) <= maxeV(i));
    end
    
    %apply the condition
    
    momXOut = momXOut(cond, :);
    momYOut = momYOut(cond, :);
    momZOut = momZOut(cond, :);
    varTof = tofOut(~cond, :);  
    cutTof = [cutTof; varTof(:)];
    tofOut = tofOut(cond, :);
    partEnergyOut = partEnergyOut(cond, :);
    hitNoOut = hitNoOut(cond, :);
    shotNoOut = shotNoOut(cond, :);
    numHitsOut = numHitsOut(cond, :);

    clear cond
    
end

%find the total momentum of each supposed coincidence
totalMomAll = ((sum(momXOut.*repmat(incl, size(momXOut,1),1), 2)).^2+...
    (sum(momYOut.*repmat(incl, size(momXOut,1),1), 2)).^2+...
    (sum(momZOut.*repmat(incl, size(momXOut,1),1), 2)).^2).^(1/2);

%find the energy of the coincidence
KER = sum(partEnergyOut, 2);

%looking at t1 v t2
if conDeltat
    switch numPart
        case 1
            cond = true(size(momXOut(:, 1)));
            
        case 2
            cond = conMomComp2(momZOut(:, 1), momZOut(:, 2), momXOut(:, 1),...
                momXOut(:, 2), momYOut(:, 1), momYOut(:, 2),...
                deltat, deltax, deltay);

        case 3
            cond = conMomComp3(momXOut(:, 1), momXOut(:, 2), momXOut(:, 3), deltax);
            cond = cond & conMomComp3(momYOut(:, 1), momYOut(:, 2), momYOut(:, 3), deltay);
            cond = cond & conMomComp3(momZOut(:, 1), momZOut(:, 2), momZOut(:, 3), deltat);
            
        case 4
            cond = conMomComp4(momZOut(:, 1), momZOut(:, 2), momZOut(:, 3), momZOut(:, 4), deltat);
            cond = cond & conMomComp4(momYOut(:, 1), momYOut(:, 2), momYOut(:, 3), momYOut(:, 4), deltay);
            cond = cond & conMomComp4(momXOut(:, 1), momXOut(:, 2), momXOut(:, 3), momXOut(:, 4), deltax);

    end
else
    cond = true(size(momXOut(:,1)));
end

if plTofPlot
    switch numPart
        case 1
        case 2
            figure()
            plot(tofOut(cond, 1), tofOut(cond, 2), 'bo',...
                tofOut(~cond, 1), tofOut(~cond, 2), 'go', 'MarkerSize', 2)
            xlabel(['tof (ns) mass ' num2str(mass(1)) ' charge ' num2str(charge(1))])
            ylabel(['tof (ns) mass ' num2str(mass(2)) ' charge ' num2str(charge(2))])
            
        case 3
            figure()
            plot(tofOut(cond, 1)+tofOut(cond, 2), tofOut(cond, 3), 'bo',...
                tofOut(~cond, 1)+tofOut(~cond, 2), tofOut(~cond, 3), 'go', 'MarkerSize', 2)
            xlabel(['tof (ns) mass ' num2str(mass(1)) ' charge ' num2str(charge(1)) ' + mass ' num2str(mass(2)) ' charge ' num2str(charge(2))])
            ylabel(['tof (ns) mass ' num2str(mass(3)) ' charge ' num2str(charge(3))])

        case 4

            figure()
            plot(tofOut(cond, 1)+tofOut(cond, 3), tofOut(cond, 2)+tofOut(cond, 4), 'bo',...
                tofOut(~cond, 1)+tofOut(~cond, 3), tofOut(~cond, 2)+tofOut(~cond, 4), 'go', 'MarkerSize', 2)
            xlabel(['tof (ns) mass ' num2str(mass(1)) ' charge ' num2str(charge(1)) ' + mass ' num2str(mass(3)) ' charge ' num2str(charge(3))])
            ylabel(['tof (ns) mass ' num2str(mass(2)) ' charge ' num2str(charge(2)) ' + mass ' num2str(mass(4)) ' charge ' num2str(charge(4))])

    end
end

%apply the condition

momXOut = momXOut(cond, :);
momYOut = momYOut(cond, :);
momZOut = momZOut(cond, :);
varTof = tofOut(~cond, :);
tofOut = tofOut(cond, :);
cutTof = [cutTof; varTof(:)];
hitNoOut = hitNoOut(cond, :);
shotNoOut = shotNoOut(cond, :);
numHitsOut = numHitsOut(cond, :);
KER = KER(cond);
partEnergyOut = partEnergyOut(cond, :);

%determine how many coincidences there are
stats(3) = size(momXOut,1);
stats(4) = length(unique(shotNoOut(:, 1)));

%remove the first 0 of cutTof
cutTof(1) = [];

%----------------------------------plots----------------------------------%

uniqueShots = length(unique(shotNoOut(:, 1))) == length(shotNoOut(:, 1));

%plot histograms of the tof for each particle
if plTofHist && numCutPlotBins > 0
    figure()
    
    hist(tofOut, numCutPlotBins)
    xlabel('tof (ns) of each species')
    title(['num points ' num2str(numel(tofOut)/numel(incl))])
    
    leg = {};
    for i = 1:length(mass)
        var = {['mass ' num2str(mass(i)) ' charge ' num2str(charge(i))]};
        leg = [leg; var];
    end
    legend(leg)
    
    figure()
    hist(tofOut(:), numCutPlotBins)
    xlabel('tof (ns) of included species')
    title(['num points ' num2str(numel(tofOut))])
   
    figure()
    hist(cutTof, numCutPlotBins)
    xlabel('tof (ns) of cut species')
    title(['num points ' num2str(numel(cutTof))])
end

%histogram the energy of each particle
if plEnergyPartHist
    figure()
    leg = {};
    for i = 1:length(mass)
        var = {['mass ' num2str(mass(i)) ' charge ' num2str(charge(i))]};
        leg = [leg; var];
    end
    
    hist(partEnergyOut, numCutPlotBins)
    xlabel('energy of each species (eV)')
    title(['num points ' num2str(size(partEnergyOut, 1))])
    legend(leg)
end

%histogram the total energy of the coincidence
if plEnergyHist
    figure()
    hist(KER, numCutPlotBins)
    xlabel('total KER (eV)')
    title(['num points ' num2str(numel(KER))])
end

%plot histograms of the sum of momentum in each direction
if plMomDir

    [vals, bins] = hist((sum(momXOut.*repmat(incl, size(momXOut,1),1), 2)), numCutPlotBins);
    figure()
    plot(bins, vals)
    xlabel('x-momentum')
    
    [vals, bins] = hist((sum(momYOut.*repmat(incl, size(momXOut,1),1), 2)), numCutPlotBins);
    figure()
    plot(bins, vals)
    xlabel('y-momentum')
    
    [vals, bins] = hist((sum(momZOut.*repmat(incl, size(momXOut,1),1), 2)), numCutPlotBins);
    figure()
    plot(bins, vals)
    xlabel('z-momentum')
end

%-------------------------build output structure--------------------------%

output = struct('hitNoOut', hitNoOut, 'shotNoOut', shotNoOut, 'numHitsOut',...
    numHitsOut, 'momXOut', momXOut, 'momYOut', momYOut, 'momZOut', momZOut,...
    'cutTof', cutTof, 'tofOut', tofOut, 'mass', mass, 'charge', charge,...
    'stats', stats, 'totalMomAll', totalMomAll, 'uniqueShots',...
    uniqueShots, 'partEnergyOut', partEnergyOut, 'KER', KER, 'incl', incl);



function [momXOut, momYOut, momZOut, partEnergyOut, tofOut,...
          hitNoOut,shotNoOut, numHitsOut, cutTof] = ...
          FirstCut(numPart, numExtraHits, colms, momX, momY, momZ, ...
                   partEnergy, tof, hitNo, shotNo, numHits)
%first cut: require that the total number of hits equals the total number
%of particles in the desired coincidence
expectLength = 0;

cutTof = 0;
momXOut =[];
momYOut = [];
momZOut = [];
partEnergyOut = [];
tofOut = [];
hitNoOut = [];
shotNoOut = [];
numHitsOut = [];
condd = false(size(tof));
        
for j = numPart:(numExtraHits+numPart)

    last_indeces = find((hitNo == j) & (numHits == j));
    combs = nchoosek(1:j, numPart);

    for i = 1:size(combs, 1)

        momXOut_i = [];
        momYOut_i = [];
        momZOut_i = [];
        partEnergy_i = [];
        tofOut_i = [];
        hitNoOut_i = [];
        shotNoOut_i = [];
        numHitsOut_i = [];

        for k = 1:numPart

            particle_indeces = last_indeces - combs(i, numPart - k + 1) + 1;
            condd(particle_indeces) = true;
            momXOut_i = [momXOut_i, momX(particle_indeces, colms(k))];
            momYOut_i = [momYOut_i, momY(particle_indeces, colms(k))];
            momZOut_i = [momZOut_i, momZ(particle_indeces, colms(k))];
            partEnergy_i = [partEnergy_i, partEnergy(particle_indeces, colms(k))];
            tofOut_i = [tofOut_i, tof(particle_indeces)];
            hitNoOut_i = [hitNoOut_i, hitNo(particle_indeces)];
            shotNoOut_i = [shotNoOut_i, shotNo(particle_indeces)];
            numHitsOut_i = [numHitsOut_i, numHits(particle_indeces)];

        end
        momXOut = [momXOut; momXOut_i];
        momYOut = [momYOut; momYOut_i];
        momZOut = [momZOut; momZOut_i];
        partEnergyOut = [partEnergyOut; partEnergy_i];
        tofOut = [tofOut; tofOut_i];
        hitNoOut = [hitNoOut; hitNoOut_i];
        shotNoOut = [shotNoOut; shotNoOut_i];
        numHitsOut = [numHitsOut; numHitsOut_i];
    end
end

cutTof = [cutTof; tof(~condd)];