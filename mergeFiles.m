function mergeFiles(files, filename, numParts, multiplicity, deleteFiles)

clear output

hitNoOut = [];
shotNoOut = [];
numHitsOut = [];
momXOut = [];
momYOut = [];
momZOut = [];
tofOut = [];
mass = [];
charge = [];
maxeV = [];
partEnergy = [];
KER = [];
angle = [];
heavyAngle = [];
incl = [];
mom = [];
cutTof = [];
totalMom = [];
totalMomCut = [];
totalMomAll = [];
uniqueShots = [];
angleLog = [];
numExtraHits = [];
coneV = [];
conTotMom = [];
maxMom = [];
minMom = [];
conHeavyPartAngle = [];
maxHeavyPartAngle = [];
stats = [];
pxc = [];
pyc = [];
pzc = [];
V1 = [];
VM = [];
s = [];
A = [];
L = [];
C = [];
conTotEnergy = [];
minKER = [];
maxKER = [];
aa = [];
bb = [];
cc = [];
methods = [];
t0 = [];
x0 = [];

for ii=1:length(files)

    load(files{ii})
    
    hitNo = output.hitNoOut;
    shotNo = output.shotNoOut;
    numHits = output.numHitsOut;
    momX = output.momXOut;
    momY = output.momYOut;
    momZ = output.momZOut;
    tof = output.tofOut;
    mass1 = output.mass;
    charge1 = output.charge;
    partEnergy1 = output.partEnergy;
    KER1 = output.KER;
    angle1 = output.angle;
    heavyAngle1 = output.heavyAngle;
    incl1 = output.incl;
    mom1 = output.mom;
    maxeV1 = output.maxeV;
    cutTof1 = output.cutTof;
    totalMom1 = output.totalMom;
    totalMomCut1 = output.totalMomCut;
    totalMomAll1 = output.totalMomAll;
    uniqueShots1 = output.uniqueShots;
    angleLog1 = output.angleLog;
    numExtraHits1 = output.numExtraHits;
    coneV1 = output.coneV;
    conTotMom1 = output.conTotMom;
    maxMom1 = output.maxMom;
    minMom1 = output.minMom;
    conHeavyPartAngle1 = output.conHeavyPartAngle;
    maxHeavyPartAngle1 = output.maxHeavyPartAngle;
    stats1 = output.stats;
    pxc1 = output.pxc;
    pyc1 = output.pyc;
    pzc1 = output.pzc;
    V11 = output.V1;
    VM1 = output.VM;
    s1 = output.s;
    A1 = output.A;
    L1 = output.L;
    C1 = output.C;
    conTotEnergy1 = output.conTotEnergy;
    minKER1 = output.minKER;
    maxKER1 = output.maxKER;
    aa1 = output.aa;
    bb1 = output.bb;
    cc1 = output.cc;
    methods1 = output.methods;
    t01 = output.t0;
    x01 = output.x0;
        
    [~, R] = sort(mass1./charge1);
    
    hitNo = hitNo(:,R);
    shotNo = shotNo(:,R);
    numHits = numHits(:,R);
    momX = momX(:,R);
    momY = momY(:,R);
    momZ = momZ(:,R);
    tof = tof(:,R);
    mass1 = mass1(:,R);
    charge1 = charge1(:,R);
    partEnergy1 = partEnergy1(:,R);
    incl1 = incl1(:,R);
    mom1 = mom1(:,R);
    maxeV1 = maxeV1(:,R);
    
    [~, R] = sort(incl1, 'descend');
    
    hitNo = hitNo(:,R);
    shotNo = shotNo(:,R);
    numHits = numHits(:,R);
    momX = momX(:,R);
    momY = momY(:,R);
    momZ = momZ(:,R);
    tof = tof(:,R);
    mass1 = mass1(:,R);
    charge1 = charge1(:,R);
    partEnergy1 = partEnergy1(:,R);
    incl1 = incl1(:,R);
    mom1 = mom1(:,R);
    maxeV1 = maxeV1(:,R);
    
    [WW,LL] = size(tof);
    
    if LL == numParts+2
    
        hitNoOut = [hitNoOut; hitNo];
        shotNoOut = [shotNoOut; shotNo];
        numHitsOut = [numHitsOut; numHits];
        momXOut = [momXOut; momX];
        momYOut = [momYOut; momY];
        momZOut = [momZOut; momZ];
        tofOut = [tofOut; tof];
        mass = [mass; mass1];
        charge = [charge; charge1];
        partEnergy = [partEnergy; partEnergy1];
        KER = [KER; KER1];
        angle = [angle; angle1];
        heavyAngle = [heavyAngle; heavyAngle1];
        incl = [incl; incl1];
        mom = [mom; mom1];
        maxeV = [maxeV; maxeV1];
        cutTof = [cutTof; cutTof1];
        totalMom = [totalMom; totalMom1];
        totalMomCut = [totalMomCut; totalMomCut1];
        totalMomAll = [totalMomAll; totalMomAll1];
        uniqueShots = [uniqueShots; uniqueShots1];
        angleLog = {angleLog, angleLog1};
        numExtraHits = [numExtraHits; numExtraHits1];
        coneV = [coneV; coneV1];
        conTotMom = [conTotMom; conTotMom1];
        maxMom = [maxMom; maxMom1];
        minMom = [minMom; minMom1];
        conHeavyPartAngle = [conHeavyPartAngle; conHeavyPartAngle1];
        maxHeavyPartAngle = [maxHeavyPartAngle; maxHeavyPartAngle1];
        stats = {stats; stats1};
        pxc = [pxc; pxc1];
        pyc = [pyc; pyc1];
        pzc = [pzc; pzc1];
        V1 = [V1; V11];
        VM = [VM; VM1];
        s = [s; s1];
        A = [A; A1];
        L = [L; L1];
        C = [C; C1];
        conTotEnergy = [conTotEnergy; conTotEnergy1];
        minKER = [minKER; minKER1];
        maxKER = [maxKER; maxKER1];
        aa = [aa; aa1];
        bb = [bb; bb1];
        cc = [cc; cc1];
        methods = [methods; methods1];
        t0 = [t0; t01];
        x0 = [x0; x01];

    end
    
    if LL == numParts+1
        
        if multiplicity == 1
    
            hitNoOut = [hitNoOut; hitNo];
            shotNoOut = [shotNoOut; shotNo];
            numHitsOut = [numHitsOut; numHits];
            momXOut = [momXOut; momX];
            momYOut = [momYOut; momY];
            momZOut = [momZOut; momZ];
            tofOut = [tofOut; tof];
            mass = [mass; mass1];
            charge = [charge; charge1];
            partEnergy = [partEnergy; partEnergy1];
            KER = [KER; KER1];
            angle = [angle; angle1];
            heavyAngle = [heavyAngle; heavyAngle1];
            incl = [incl; incl1];
            mom = [mom; mom1];
            maxeV = [maxeV; maxeV1];
            cutTof = [cutTof; cutTof1];
            totalMom = [totalMom; totalMom1];
            totalMomCut = [totalMomCut; totalMomCut1];
            totalMomAll = [totalMomAll; totalMomAll1];
            uniqueShots = [uniqueShots; uniqueShots1];
            angleLog = {angleLog, angleLog1};
            numExtraHits = [numExtraHits; numExtraHits1];
            coneV = [coneV; coneV1];
            conTotMom = [conTotMom; conTotMom1];
            maxMom = [maxMom; maxMom1];
            minMom = [minMom; minMom1];
            conHeavyPartAngle = [conHeavyPartAngle; conHeavyPartAngle1];
            maxHeavyPartAngle = [maxHeavyPartAngle; maxHeavyPartAngle1];
            stats = {stats; stats1};
            pxc = [pxc; pxc1];
            pyc = [pyc; pyc1];
            pzc = [pzc; pzc1];
            V1 = [V1; V11];
            VM = [VM; VM1];
            s = [s; s1];
            A = [A; A1];
            L = [L; L1];
            C = [C; C1];
            conTotEnergy = [conTotEnergy; conTotEnergy1];
            minKER = [minKER; minKER1];
            maxKER = [maxKER; maxKER1];
            aa = [aa; aa1];
            bb = [bb; bb1];
            cc = [cc; cc1];
            methods = [methods; methods1];
            t0 = [t0; t01];
            x0 = [x0; x01];
        
        end
        
        if multiplicity == 2
            hitNoOut = [hitNoOut; [hitNo,zeros(WW,1)]];
            shotNoOut = [shotNoOut; [shotNo,zeros(WW,1)]];
            numHitsOut = [numHitsOut; [numHits,zeros(WW,1)]];
            momXOut = [momXOut; [momX,zeros(WW,1)]];
            momYOut = [momYOut; [momY,zeros(WW,1)]];
            momZOut = [momZOut; [momZ,zeros(WW,1)]];
            tofOut = [tofOut; [tof,zeros(WW,1)]];
            mass = [mass; [mass1,0]];
            charge = [charge; [charge1,0]];
            partEnergy = [partEnergy; [partEnergy1,zeros(WW,1)]];
            KER = [KER; KER1];
            angle = [angle; angle1];
            heavyAngle = [heavyAngle; heavyAngle1];
            incl = [incl; [incl1,0]];
            mom = [mom; [mom1,zeros(WW,1)]];
            maxeV = [maxeV; [maxeV1,0]];
            cutTof = [cutTof; cutTof1];
            totalMom = [totalMom; totalMom1];
            totalMomCut = [totalMomCut; totalMomCut1];
            totalMomAll = [totalMomAll; totalMomAll1];
            uniqueShots = [uniqueShots; uniqueShots1];
            angleLog = {angleLog, angleLog1};
            numExtraHits = [numExtraHits; numExtraHits1];
            coneV = [coneV; coneV1];
            conTotMom = [conTotMom; conTotMom1];
            maxMom = [maxMom; maxMom1];
            minMom = [minMom; minMom1];
            conHeavyPartAngle = [conHeavyPartAngle; conHeavyPartAngle1];
            maxHeavyPartAngle = [maxHeavyPartAngle; maxHeavyPartAngle1];
            stats = {stats; stats1};
            pxc = [pxc; pxc1];
            pyc = [pyc; pyc1];
            pzc = [pzc; pzc1];
            V1 = [V1; V11];
            VM = [VM; VM1];
            s = [s; s1];
            A = [A; A1];
            L = [L; L1];
            C = [C; C1];
            conTotEnergy = [conTotEnergy; conTotEnergy1];
            minKER = [minKER; minKER1];
            maxKER = [maxKER; maxKER1];
            aa = [aa; aa1];
            bb = [bb; bb1];
            cc = [cc; cc1];
            methods = [methods; methods1];
            t0 = [t0; t01];
            x0 = [x0; x01];
        end
    
    elseif LL == numParts
        
        if multiplicity == 1
            hitNoOut = [hitNoOut; [hitNo,zeros(WW,1)]];
            shotNoOut = [shotNoOut; [shotNo,zeros(WW,1)]];
            numHitsOut = [numHitsOut; [numHits,zeros(WW,1)]];
            momXOut = [momXOut; [momX,zeros(WW,1)]];
            momYOut = [momYOut; [momY,zeros(WW,1)]];
            momZOut = [momZOut; [momZ,zeros(WW,1)]];
            tofOut = [tofOut; [tof,zeros(WW,1)]];
            mass = [mass; [mass1,0]];
            charge = [charge; [charge1,0]];
            partEnergy = [partEnergy; [partEnergy1,zeros(WW,1)]];
            KER = [KER; KER1];
            angle = [angle; angle1];
            heavyAngle = [heavyAngle; heavyAngle1];
            incl = [incl; [incl1,0]];
            mom = [mom; [mom1,zeros(WW,1)]];
            maxeV = [maxeV; [maxeV1,0]];
            cutTof = [cutTof; cutTof1];
            totalMom = [totalMom; totalMom1];
            totalMomCut = [totalMomCut; totalMomCut1];
            totalMomAll = [totalMomAll; totalMomAll1];
            uniqueShots = [uniqueShots; uniqueShots1];
            angleLog = {angleLog, angleLog1};
            numExtraHits = [numExtraHits; numExtraHits1];
            coneV = [coneV; coneV1];
            conTotMom = [conTotMom; conTotMom1];
            maxMom = [maxMom; maxMom1];
            minMom = [minMom; minMom1];
            conHeavyPartAngle = [conHeavyPartAngle; conHeavyPartAngle1];
            maxHeavyPartAngle = [maxHeavyPartAngle; maxHeavyPartAngle1];
            stats = {stats; stats1};
            pxc = [pxc; pxc1];
            pyc = [pyc; pyc1];
            pzc = [pzc; pzc1];
            V1 = [V1; V11];
            VM = [VM; VM1];
            s = [s; s1];
            A = [A; A1];
            L = [L; L1];
            C = [C; C1];
            conTotEnergy = [conTotEnergy; conTotEnergy1];
            minKER = [minKER; minKER1];
            maxKER = [maxKER; maxKER1];
            aa = [aa; aa1];
            bb = [bb; bb1];
            cc = [cc; cc1];
            methods = [methods; methods1];
            t0 = [t0; t01];
            x0 = [x0; x01];
        end
           
        if multiplicity == 2
            hitNoOut = [hitNoOut; [hitNo,zeros(WW,1),zeros(WW,1)]];
            shotNoOut = [shotNoOut; [shotNo,zeros(WW,1),zeros(WW,1)]];
            numHitsOut = [numHitsOut; [numHits,zeros(WW,1),zeros(WW,1)]];
            momXOut = [momXOut; [momX,zeros(WW,1),zeros(WW,1)]];
            momYOut = [momYOut; [momY,zeros(WW,1),zeros(WW,1)]];
            momZOut = [momZOut; [momZ,zeros(WW,1),zeros(WW,1)]];
            tofOut = [tofOut; [tof,zeros(WW,1),zeros(WW,1)]];
            mass = [mass; [mass1,0,0]];
            charge = [charge; [charge1,0,0]];
            partEnergy = [partEnergy; [partEnergy1,zeros(WW,1),zeros(WW,1)]];
            KER = [KER; KER1];
            angle = [angle; angle1];
            heavyAngle = [heavyAngle; heavyAngle1];
            incl = [incl; [incl1,0,0]];
            mom = [mom; [mom1,zeros(WW,1),zeros(WW,1)]];
            maxeV = [maxeV; [maxeV1,0,0]];
            cutTof = [cutTof; cutTof1];
            totalMom = [totalMom; totalMom1];
            totalMomCut = [totalMomCut; totalMomCut1];
            totalMomAll = [totalMomAll; totalMomAll1];
            uniqueShots = [uniqueShots; uniqueShots1];
            angleLog = {angleLog, angleLog1};
            numExtraHits = [numExtraHits; numExtraHits1];
            coneV = [coneV; coneV1];
            conTotMom = [conTotMom; conTotMom1];
            maxMom = [maxMom; maxMom1];
            minMom = [minMom; minMom1];
            conHeavyPartAngle = [conHeavyPartAngle; conHeavyPartAngle1];
            maxHeavyPartAngle = [maxHeavyPartAngle; maxHeavyPartAngle1];
            stats = {stats; stats1};
            pxc = [pxc; pxc1];
            pyc = [pyc; pyc1];
            pzc = [pzc; pzc1];
            V1 = [V1; V11];
            VM = [VM; VM1];
            s = [s; s1];
            A = [A; A1];
            L = [L; L1];
            C = [C; C1];
            conTotEnergy = [conTotEnergy; conTotEnergy1];
            minKER = [minKER; minKER1];
            maxKER = [maxKER; maxKER1];
            aa = [aa; aa1];
            bb = [bb; bb1];
            cc = [cc; cc1];
            methods = [methods; methods1];
            t0 = [t0; t01];
            x0 = [x0; x01];
        end
        
    end

end

[~,Aa] = unique(shotNoOut(:,1));

hitNoOut = hitNoOut(Aa,:);
shotNoOut = shotNoOut(Aa,:);
numHitsOut = numHitsOut(Aa,:);
momXOut = momXOut(Aa,:);
momYOut = momYOut(Aa,:);
momZOut = momZOut(Aa,:);
tofOut = tofOut(Aa,:);
partEnergy = partEnergy(Aa,:);
KER = KER(Aa,:);
heavyAngle = heavyAngle(Aa,:);
angle = angle(Aa,:);
mom = mom(Aa,:);
totalMom = totalMom(Aa,:);

clear output

output.hitNoOut = hitNoOut;
output.shotNoOut = shotNoOut;
output.numHitsOut = numHitsOut;
output.hitNoOut = hitNoOut;
output.shotNoOut = shotNoOut;
output.numHitsOut = numHitsOut;
output.momXOut = momXOut;
output.momYOut = momYOut;
output.momZOut = momZOut;
output.tofOut = tofOut;
output.mass = mass;
output.charge = charge;
output.partEnergy = partEnergy;
output.KER = KER;
output.angle = angle;
output.heavyAngle = heavyAngle;
output.incl = incl;
output.mom = mom;
output.maxeV = maxeV;
output.cutTof = cutTof;
output.totalMom = totalMom;
output.totalMomCut = totalMomCut;
output.totalMomAll = totalMomAll;
output.uniqueShots = uniqueShots;
output.angleLog = angleLog;
output.numExtraHits = numExtraHits;
output.coneV = coneV;
output.conTotMom = conTotMom;
output.maxMom = maxMom;
output.minMom = minMom;
output.conHeavyPartAngle = conHeavyPartAngle;
output.maxHeavyPartAngle = maxHeavyPartAngle;
output.stats = stats;
output.pxc = pxc;
output.pyc = pyc;
output.pzc = pzc;
output.V1 = V1;
output.VM = VM;
output.s = s;
output.A = A;
output.L = L;
output.C = C;
output.conTotEnergy = conTotEnergy;
output.minKER = minKER;
output.maxKER = maxKER;
output.aa = aa;
output.bb = bb;
output.cc = cc;
output.methods = methods;
output.t0 = t0;
output.x0 = x0;

[path, ~] = fileparts(files{1});
filename = fullfile(path, filename);
save(filename, 'output', '-v7.3')

if deleteFiles
    for ii = 1:max(size(files))
        delete(files{ii})
    end
end
