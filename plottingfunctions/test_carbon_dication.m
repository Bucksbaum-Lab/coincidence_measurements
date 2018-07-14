histx = linspace(0,28,29);

figure()
[hist266,X266] = hist(output_266_14.partEnergyOut(output_266_14.KER<18,2),histx);
[hist1300,X1300] = hist(output_1300_14.partEnergyOut(output_1300_14.KER<18,2),histx);
[hist800,X800] = hist(output_800_14.partEnergyOut(output_800_14.KER<18,2),histx);
plot(X266,hist266/sum(hist266),X800,hist800/sum(hist800),X1300,hist1300/sum(hist1300))
legend(['266, mean eV ' num2str(mean(output_266_14.partEnergyOut(output_266_14.KER<18,2)))],...
    ['800, mean eV ' num2str(mean(output_800_14.partEnergyOut(output_800_14.KER<18,2)))],...
    ['1300, mean eV ' num2str(mean(output_1300_14.partEnergyOut(output_1300_14.KER<18,2)))])
title('CH2/C')

figure()
[hist266,X266] = hist(output_266_1.partEnergyOut(output_266_1.KER<18,2),histx);
[hist1300,X1300] = hist(output_1300_1.partEnergyOut(output_1300_1.KER<18,2),histx);
[hist800,X800] = hist(output_800_1.partEnergyOut(output_800_1.KER<18,2),histx);
plot(X266,hist266/sum(hist266),X800,hist800/sum(hist800),X1300,hist1300/sum(hist1300))
legend(['266, mean eV ' num2str(mean(output_266_1.partEnergyOut(output_266_1.KER<18,2)))],...
    ['800, mean eV ' num2str(mean(output_800_1.partEnergyOut(output_800_1.KER<18,2)))],...
    ['1300, mean eV ' num2str(mean(output_1300_1.partEnergyOut(output_1300_1.KER<18,2)))])
title('C2H/H')

figure()
[hist266,X266] = hist(output_266_13.partEnergyOut(output_266_13.KER<18,2),histx);
[hist1300,X1300] = hist(output_1300_13.partEnergyOut(output_1300_13.KER<18,2),histx);
[hist800,X800] = hist(output_800_13.partEnergyOut(output_800_13.KER<18,2),histx);
plot(X266,hist266/sum(hist266),X800,hist800/sum(hist800),X1300,hist1300/sum(hist1300))
legend(['266, mean eV ' num2str(mean(output_266_13.partEnergyOut(output_266_13.KER<18,2)))],...
    ['800, mean eV ' num2str(mean(output_800_13.partEnergyOut(output_800_13.KER<18,2)))],...
    ['1300, mean eV ' num2str(mean(output_1300_13.partEnergyOut(output_1300_13.KER<18,2)))])
title('CH/CH')