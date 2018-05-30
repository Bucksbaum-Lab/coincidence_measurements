function condition = ApplyExperimentType(condition, shutterChoice, intensityChoice,...
    paramChoice, polarChoice, delayChoice, shutterStatus, intensityStatus, paramStatus,...
    polarizationStatus, delayStatus, polarInfo, delayInfo, paramInfo)

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

for nn = 1:width(paramInfo(2,:))
    if paramChoice == nn+1
        condition = condition & (paramStatus == str2double(cell2mat(table2array(paramInfo(1,nn)))));
    end
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