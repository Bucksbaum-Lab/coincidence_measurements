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

%---------------------Last Edited on February 4, 2014---------------------%
%------Last Edited by and Original Author - Chelsea Liekhus-Schmaltz------% 
%}

function [output] = cutandplot(mass, charge, colms, momX, momY, momZ, partEnergy,...
    tof, hitNo, shotNo, numHits, numExtraHits, coneV, maxeV,...
    conHeavyPartAngle, maxHeavyPartAngle, plTofPlot,...
    plTofHist, plAngleHist, plEnergyPartHist, plEnergyHist, plMomDir,...
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
cutTof = 0;
stats = zeros(numPart+4, 1);
stats(1) = length(tof);
stats(2) = length(unique(shotNo));

%record how many of each particle there is before any cuts
for i = 5:numPart+4
    stats(i) = sum(sum(~isnan(partEnergy)));
end

%first cut: require that the total number of hits equals the total number
%of particles in the desired coincidence
if numExtraHits == 0

    %generate the condition
    cond = numHits == numPart;
    
    %apply the cut
    momX = momX(cond, :);
    momY = momY(cond, :);
    momZ = momZ(cond, :);
    hitNo = hitNo(cond);
    shotNo = shotNo(cond);
    numHits = numHits(cond);
    cutTof = [cutTof; tof(~cond)];
    tof = tof(cond);
    partEnergy = partEnergy(cond,:);

    %initialize the momentum vectors
    leng = sum(cond)/numPart;
    momXOut = zeros(leng, numPart);
    momYOut = zeros(leng, numPart);
    momZOut = zeros(leng, numPart);
    partEnergyOut = zeros(leng, numPart);
    tofOut = zeros(leng, numPart);
    hitNoOut = zeros(leng, numPart);
    shotNoOut = zeros(leng, numPart);
    numHitsOut = zeros(leng, numPart);

    %associate each hit with a mass
    for i = 1:numPart
        momXOut(:, i) = momX(hitNo == i, colms(i));
        momYOut(:, i) = momY(hitNo == i, colms(i));
        momZOut(:, i) = momZ(hitNo == i, colms(i));
        partEnergyOut(:, i) = partEnergy(hitNo == i, colms(i));
        tofOut(:, i) = tof(hitNo == i);
        hitNoOut(:, i) = hitNo(hitNo == i);
        shotNoOut(:, i) = shotNo(hitNo == i);
        numHitsOut(:, i) = numHits(hitNo == i);
    end

%if this cut is not applied then cycle through possible combinations that
%result in the necessary number of hits
else
    switch numPart
        %for one particle, the logic is straight forward, just take the
        %momentum from the appropriate column from the momentum matrix
        %generted in the preparation
        case 1
            cond = numHits <= numExtraHits+1;
            momXOut = momX(cond, colms(1));
            momYOut = momY(cond, colms(1));
            momZOut = momZ(cond, colms(1));
            partEnergyOut = partEnergy(cond, colms(1));
            tofOut = tof(cond);
            hitNoOut = hitNo(cond);
            shotNoOut = shotNo(cond);
            numHitsOut = numHits(cond);

            condd = cond;
            
            clear cond
            
        %for the two particle case, generate ordered pairs of the hits and
        %then associate these with the correct mass
        case 2
            expectLength = 0;
            
            momXOut = [];
            momYOut = [];
            momZOut = [];
            partEnergyOut = [];
            tofOut = [];
            hitNoOut = [];
            shotNoOut = [];
            numHitsOut = [];
            
            for j = 2:numExtraHits+2
                for i = 1:j-1

                    cond2 = hitNo == j;
                    cond1 = circshift(cond2, -i);

                    momXOut = [momXOut;[momX(cond1, colms(1)), momX(cond2, colms(2))]];
                    momYOut = [momYOut;[momY(cond1, colms(1)), momY(cond2, colms(2))]];
                    momZOut = [momZOut;[momZ(cond1, colms(1)), momZ(cond2, colms(2))]];
                    partEnergyOut = [partEnergyOut;[partEnergy(cond1, colms(1)), partEnergy(cond2, colms(2))]];
                    tofOut = [tofOut;[tof(cond1), tof(cond2)]];
                    hitNoOut = [hitNoOut;[hitNo(cond1), hitNo(cond2)]];
                    shotNoOut = [shotNoOut;[shotNo(cond1), shotNo(cond2)]];
                    numHitsOut = [numHitsOut;[numHits(cond1), numHits(cond2)]];

                    %record all of the events that are categorized
                    condd = cond1 | cond2;

                end
            end
            
            clear cond1
            clear cond2    
    end
    
    if isempty(condd)
        condd = false(size(tof));
    end
    
    cutTof = [cutTof; tof(~condd)];
end

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

%plot histograms of the angles
if plAngleHist && numCutPlotBins > 0
    switch numPart

        case 2
            figure()
            hist(heavyAngle, numCutPlotBins);
            xlabel('angle (radians)')
            title(['num points ' num2str(numel(heavyAngle))])
            
        case 3
            figure()
            hist(angle(:, 1), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(3)) ' charge '...
                num2str(charge(3)) ' and mass ' num2str(mass(1)) ' charge ' num2str(charge(1))])
            title(['num points ' num2str(size(angle,1))])
            
            figure()
            hist(angle(:, 2), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(3)) ' charge ' num2str(charge(3))...
                ' and mass ' num2str(mass(2)) ' charge ' num2str(charge(2))])
            title(['num points ' num2str(size(angle,1))])
            
            figure()
            hist(angle(:, 3), numCutPlotBins)
            xlabel(['angle between center of mass (' num2str(mass(2)) ' charge '...
                num2str(charge(2)) ' - mass ' num2str(mass(3)) ' charge ' num2str(charge(3))...
                ') and mass ' num2str(mass(1)) ' charge ' num2str(charge(1))])
            title(['num points ' num2str(size(angle,1))])
            
        case 4
            figure()
            hist(angle(:, 1), numCutPlotBins)
            xlabel(['angle between center of mass (' num2str(mass(3)) ' charge '...
                num2str(charge(3)) ' - mass ' num2str(mass(3)) 'charge ' num2str(charge(4))...
                ') and mass ' num2str(mass(1)) ' charge ' num2str(charge(1))])
            title(['num points ' num2str(size(angle,1))])
            
            figure()
            hist(angle(:, 2), numCutPlotBins)
            xlabel(['angle between center of mass (' num2str(mass(3)) ' charge '...
                num2str(charge(3)) ' - mass ' num2str(mass(3)) 'charge ' num2str(charge(4))...
                ') and mass ' num2str(mass(2)) ' charge ' num2str(charge(2))])
            title(['num points ' num2str(size(angle,1))])
            
            figure()
            hist(angle(:, 3), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(1)) ' charge '...
                num2str(charge(1)) ' and mass ' num2str(mass(2)) ' charge ' num2str(charge(2))])
            title(['num points ' num2str(size(angle,1))])
        
            figure()
            hist(angle(:, 4), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(1)) ' charge ' num2str(charge(1))...
                ' and mass ' num2str(mass(3)) ' charge ' num2str(charge(3))])
            title(['num points ' num2str(size(angle,1))])
        
            figure()
            hist(angle(:, 5), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(1)) ' charge ' num2str(charge(1))...
                ' and mass ' num2str(mass(4)) ' charge ' num2str(charge(4))])
            title(['num points ' num2str(size(angle,1))])
        
            figure()
            hist(angle(:, 6), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(2)) ' charge ' num2str(charge(2))...
                ' and mass ' num2str(mass(3)) ' charge ' num2str(charge(3))])
            title(['num points ' num2str(size(angle,1))])
        
            figure()
            hist(angle(:, 7), numCutPlotBins)
            xlabel(['angle between mass ' num2str(mass(2)) ' charge ' num2str(charge(2))...
                ' and mass ' num2str(mass(4)) ' charge ' num2str(charge(4))])
            title(['num points ' num2str(size(angle,1))])
    end
end

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