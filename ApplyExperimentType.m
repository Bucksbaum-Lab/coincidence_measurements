function condition = ApplyExperimentType(condition, shutterChoice, intensityChoice,...
    paramChoice, polarChoice, delayChoice, shutterStatus, intensityStatus, paramStatus,...
    polarizationStatus, delayStatus, polarInfo, delayInfo)

if shutterChoice == 2
    condition = condition & ~shutterStatus;
elseif shutterChoice == 3
    condition = condition & shutterStatus;
end

if intensityChoice == 2
    condition = condition & ~intensityStatus;
elseif intensityChoice == 3
    condition = condition & intensityStatus;
end

if paramChoice == 2
    condition = condition & (paramStatus == 1);
elseif paramChoice == 3
    condition = condition & (paramStatus == 2);
elseif paramChoice == 4
    condition = condition & (paramStatus == 3);
elseif paramChoice == 5
    condition = condition & (paramStatus == 4);
end

for nn = 1:length(polarInfo(2,:))
    if polarChoice == nn+1
        condition = condition & (polarizationStatus == polarInfo(2,nn));
    end
end

for nn = 1:length(delayInfo(2,:))
    if delayChoice == nn+1
        condition = condition & (delayStatus == delayInfo(2,nn));
    end
end