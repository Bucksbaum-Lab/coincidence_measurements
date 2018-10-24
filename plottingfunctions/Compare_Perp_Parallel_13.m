histx = linspace(0,18,18*1+1);

CHeV = 15;
mom = 10;
totKER = 18;

coinNum = linspace(1,size(output_266_13.numHitsOut,1),size(output_266_13.numHitsOut,1));
coinNum = coinNum';
cond_266 = abs(sum(output_266_13.momXOut,2))<mom&...
    abs(sum(output_266_13.momYOut,2))<mom&...
    abs(sum(output_266_13.momZOut,2))<mom&...
    output_266_13.KER<totKER&...
    output_266_13.partEnergyOut(:,1)<CHeV&output_266_13.partEnergyOut(:,2)<CHeV;

%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_1300_13_perp.numHitsOut,1),size(output_1300_13_perp.numHitsOut,1));
coinNum = coinNum';
cond_1300_13_perp = abs(sum(output_1300_13_perp.momXOut,2))<mom&...
    abs(sum(output_1300_13_perp.momYOut,2))<mom&...
    abs(sum(output_1300_13_perp.momZOut,2))<mom&...
    output_1300_13_perp.KER<totKER&...
    output_1300_13_perp.partEnergyOut(:,1)<CHeV&output_1300_13_perp.partEnergyOut(:,2)<CHeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_1300_13_par.numHitsOut,1),size(output_1300_13_par.numHitsOut,1));
coinNum = coinNum';
cond_1300_13_par = abs(sum(output_1300_13_par.momXOut,2))<mom&...
    abs(sum(output_1300_13_par.momYOut,2))<mom&...
    abs(sum(output_1300_13_par.momZOut,2))<mom&...
    output_1300_13_par.KER<totKER&...
    output_1300_13_par.partEnergyOut(:,1)<CHeV&output_1300_13_par.partEnergyOut(:,2)<CHeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_800_13_perp.numHitsOut,1),size(output_800_13_perp.numHitsOut,1));
coinNum = coinNum';
cond_800_13_perp = abs(sum(output_800_13_perp.momXOut,2))<mom&...
    abs(sum(output_800_13_perp.momYOut,2))<mom&...
    abs(sum(output_800_13_perp.momZOut,2))<mom&...
    output_800_13_perp.KER<totKER&...
    output_800_13_perp.partEnergyOut(:,1)<CHeV&output_800_13_perp.partEnergyOut(:,2)<CHeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_800_13_par.numHitsOut,1),size(output_800_13_par.numHitsOut,1));
coinNum = coinNum';
cond_800_13_par = abs(sum(output_800_13_par.momXOut,2))<mom&...
    abs(sum(output_800_13_par.momYOut,2))<mom&...
    abs(sum(output_800_13_par.momZOut,2))<mom&...
    output_800_13_par.KER<totKER&...
    output_800_13_par.partEnergyOut(:,1)<CHeV&output_800_13_par.partEnergyOut(:,2)<CHeV;
%&...
%    coinNum<max(coinNum)/2;

[hist266,X266] = hist(output_266_13.KER(cond_266),histx);
[hist1300_perp,X1300] = hist(output_1300_13_perp.KER(cond_1300_13_perp),histx);
[hist800_perp,X800] = hist(output_800_13_perp.KER(cond_800_13_perp),histx);
[hist1300_par,X1300] = hist(output_1300_13_par.KER(cond_1300_13_par),histx);
[hist800_par,X800] = hist(output_800_13_par.KER(cond_800_13_par),histx);

error_266 = sqrt(hist266)/2;

figure()
errorbar(X266,hist266/2,sqrt(hist266)/2,sqrt(hist266)/2)
hold on;

error_800_perp = sqrt(hist800_perp);
error_800_par = sqrt(hist800_par);

errorbar(X800,hist800_perp,sqrt(hist800_perp),sqrt(hist800_perp))
errorbar(X800,hist800_par,sqrt(hist800_par),sqrt(hist800_par))
xlabel('KER (eV)')
legend('266','800 perp','800 par')
title('CH/CH 800')

figure()
errorbar(X266,hist266/2,sqrt(hist266)/2,sqrt(hist266)/2)
hold on;

error_1300_perp = sqrt(hist1300_perp);
error_1300_par = sqrt(hist1300_par);

errorbar(X1300,hist1300_perp,sqrt(hist1300_perp),sqrt(hist1300_perp))
errorbar(X1300,hist1300_par,sqrt(hist1300_par),sqrt(hist1300_par))
xlabel('KER (eV)')
legend('266','1300 perp','1300 par')
title('CH/CH 1300')



%%

figure()
errorbar(X800,(hist800_perp-hist266/2),sqrt(error_800_perp.^2+error_266.^2),sqrt(error_800_perp.^2+error_266.^2))
hold on
errorbar(X800,(hist800_par-hist266/2),sqrt(error_800_par.^2+error_266.^2),sqrt(error_800_par.^2+error_266.^2))
legend('800 perp', '800 par')
title('CH/CH 800')
xlabel('eV')
ylabel('change in counts from 266')


figure()
errorbar(X1300,(hist1300_perp-hist266/2),sqrt(error_1300_perp.^2+error_266.^2),sqrt(error_1300_perp.^2+error_266.^2))
hold on
errorbar(X1300,(hist1300_par-hist266/2),sqrt(error_1300_par.^2+error_266.^2),sqrt(error_1300_par.^2+error_266.^2))
legend('1300 perp', '1300 par')
title('CH/CH 1300')
xlabel('eV')
ylabel('change in counts from 266')

%%
BondMomX = output_266_13.momXOut(cond_266,1)-output_266_13.momXOut(cond_266,2);
BondMomY = output_266_13.momYOut(cond_266,1)-output_266_13.momYOut(cond_266,2);
BondMomZ = output_266_13.momZOut(cond_266,1)-output_266_13.momZOut(cond_266,2);

KER = output_266_13.KER(cond_266);

tt = '266';
sp = 'symmetric';

KERAnglePlot


BondMomX = output_800_13_perp.momXOut(cond_800_13_perp,1)-output_800_13_perp.momXOut(cond_800_13_perp,2);
BondMomY = output_800_13_perp.momYOut(cond_800_13_perp,1)-output_800_13_perp.momYOut(cond_800_13_perp,2);
BondMomZ = output_800_13_perp.momZOut(cond_800_13_perp,1)-output_800_13_perp.momZOut(cond_800_13_perp,2);

KER = output_800_13_perp.KER(cond_800_13_perp);

tt = '800 perp';
sp = 'symmetric';

KERAnglePlot


BondMomX = output_1300_13_perp.momXOut(cond_1300_13_perp,1)-output_1300_13_perp.momXOut(cond_1300_13_perp,2);
BondMomY = output_1300_13_perp.momYOut(cond_1300_13_perp,1)-output_1300_13_perp.momYOut(cond_1300_13_perp,2);
BondMomZ = output_1300_13_perp.momZOut(cond_1300_13_perp,1)-output_1300_13_perp.momZOut(cond_1300_13_perp,2);

KER = output_1300_13_perp.KER(cond_1300_13_perp);

tt = '1300 perp';
sp = 'symmetric';

KERAnglePlot


BondMomX = output_800_13_par.momXOut(cond_800_13_par,1)-output_800_13_par.momXOut(cond_800_13_par,2);
BondMomY = output_800_13_par.momYOut(cond_800_13_par,1)-output_800_13_par.momYOut(cond_800_13_par,2);
BondMomZ = output_800_13_par.momZOut(cond_800_13_par,1)-output_800_13_par.momZOut(cond_800_13_par,2);

KER = output_800_13_par.KER(cond_800_13_par);

tt = '800 par';
sp = 'symmetric';

KERAnglePlot


BondMomX = output_1300_13_par.momXOut(cond_1300_13_par,1)-output_1300_13_par.momXOut(cond_1300_13_par,2);
BondMomY = output_1300_13_par.momYOut(cond_1300_13_par,1)-output_1300_13_par.momYOut(cond_1300_13_par,2);
BondMomZ = output_1300_13_par.momZOut(cond_1300_13_par,1)-output_1300_13_par.momZOut(cond_1300_13_par,2);

KER = output_1300_13_par.KER(cond_1300_13_par);

tt = '1300 par';
sp = 'symmetric';

KERAnglePlot
