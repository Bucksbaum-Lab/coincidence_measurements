histx = linspace(0,18,19);

figure()
[hist266,X266] = hist(output_266_14.KER(output_266_14.KER<18),histx);
[hist1300,X1300] = hist(output_1300_14.KER(output_1300_14.KER<18),histx);
[hist800,X800] = hist(output_800_14.KER(output_800_14.KER<18),histx);
plot(X266,hist266/sum(hist266),'-o',X800,hist800/sum(hist800),'-o',X1300,hist1300/sum(hist1300),'-o')
legend(['266, mean eV ' num2str(mean(output_266_14.KER(output_266_14.KER<18)))],...
    ['800, mean eV ' num2str(mean(output_800_14.KER(output_800_14.KER<18)))],...
    ['1300, mean eV ' num2str(mean(output_1300_14.KER(output_1300_14.KER<18)))])
title('CH2/C')

figure()
[hist266,X266] = hist(output_266_1.KER(output_266_1.KER<18),histx);
[hist1300,X1300] = hist(output_1300_1.KER(output_1300_1.KER<18),histx);
[hist800,X800] = hist(output_800_1.KER(output_800_1.KER<18),histx);
plot(X266,hist266/sum(hist266),'-o',X800,hist800/sum(hist800),'-o',X1300,hist1300/sum(hist1300), '-o')
legend(['266, mean eV ' num2str(mean(output_266_1.KER(output_266_1.KER<18)))],...
    ['800, mean eV ' num2str(mean(output_800_1.KER(output_800_1.KER<18)))],...
    ['1300, mean eV ' num2str(mean(output_1300_1.KER(output_1300_1.KER<18)))])
title('C2H/H')

figure()
[hist266,X266] = hist(output_266_13.KER(output_266_13.KER<18),histx);
[hist1300,X1300] = hist(output_1300_13.KER(output_1300_13.KER<18),histx);
[hist800,X800] = hist(output_800_13.KER(output_800_13.KER<18),histx);
plot(X266,hist266/sum(hist266),'-o',X800,hist800/sum(hist800),'-o',X1300,hist1300/sum(hist1300), '-o')
legend(['266, mean eV ' num2str(mean(output_266_13.KER(output_266_13.KER<18)))],...
    ['800, mean eV ' num2str(mean(output_800_13.KER(output_800_13.KER<18)))],...
    ['1300, mean eV ' num2str(mean(output_1300_13.KER(output_1300_13.KER<18)))])
title('CH/CH')