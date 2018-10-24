histx = linspace(0,18,18*3+1);

C2HeV = 2;
HeV = 15;
mom = 9;
totKER = 18;

coinNum = linspace(1,size(output_266_1.numHitsOut,1),size(output_266_1.numHitsOut,1));
coinNum = coinNum';
cond_266 = abs(sum(output_266_1.momXOut,2))<mom&...
    abs(sum(output_266_1.momYOut,2))<mom&...
    abs(sum(output_266_1.momZOut,2))<mom&...
    output_266_1.KER<totKER&...
    output_266_1.partEnergyOut(:,1)<HeV&output_266_1.partEnergyOut(:,2)<C2HeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_1300_1_perp.numHitsOut,1),size(output_1300_1_perp.numHitsOut,1));
coinNum = coinNum';
cond_1300_1_perp = abs(sum(output_1300_1_perp.momXOut,2))<mom&...
    abs(sum(output_1300_1_perp.momYOut,2))<mom&...
    abs(sum(output_1300_1_perp.momZOut,2))<mom&...
    output_1300_1_perp.KER<totKER&...
    output_1300_1_perp.partEnergyOut(:,1)<HeV&output_1300_1_perp.partEnergyOut(:,2)<C2HeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_1300_1_par.numHitsOut,1),size(output_1300_1_par.numHitsOut,1));
coinNum = coinNum';
cond_1300_1_par = abs(sum(output_1300_1_par.momXOut,2))<mom&...
    abs(sum(output_1300_1_par.momYOut,2))<mom&...
    abs(sum(output_1300_1_par.momZOut,2))<mom&...
    output_1300_1_par.KER<totKER&...
    output_1300_1_par.partEnergyOut(:,1)<HeV&output_1300_1_par.partEnergyOut(:,2)<C2HeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_800_1_perp.numHitsOut,1),size(output_800_1_perp.numHitsOut,1));
coinNum = coinNum';
cond_800_1_perp = abs(sum(output_800_1_perp.momXOut,2))<mom&...
    abs(sum(output_800_1_perp.momYOut,2))<mom&...
    abs(sum(output_800_1_perp.momZOut,2))<mom&...
    output_800_1_perp.KER<totKER&...
    output_800_1_perp.partEnergyOut(:,1)<HeV&output_800_1_perp.partEnergyOut(:,2)<C2HeV;
%&...
%    coinNum<max(coinNum)/2;

coinNum = linspace(1,size(output_800_1_par.numHitsOut,1),size(output_800_1_par.numHitsOut,1));
coinNum = coinNum';
cond_800_1_par = abs(sum(output_800_1_par.momXOut,2))<mom&...
    abs(sum(output_800_1_par.momYOut,2))<mom&...
    abs(sum(output_800_1_par.momZOut,2))<mom&...
    output_800_1_par.KER<totKER&...
    output_800_1_par.partEnergyOut(:,1)<HeV&output_800_1_par.partEnergyOut(:,2)<C2HeV;
%&...
%    coinNum<max(coinNum)/2;

KER_266_1 = getCOMKER2(output_266_1.momXOut(cond_266,1),output_266_1.momXOut(cond_266,2),...
    output_266_1.momYOut(cond_266,1),output_266_1.momYOut(cond_266,2),...
    output_266_1.momZOut(cond_266,1),output_266_1.momZOut(cond_266,2),1,25);

KER_1300_perp_1 = getCOMKER2(output_1300_1_perp.momXOut(cond_1300_1_perp,1),output_1300_1_perp.momXOut(cond_1300_1_perp,2),...
    output_1300_1_perp.momYOut(cond_1300_1_perp,1),output_1300_1_perp.momYOut(cond_1300_1_perp,2),...
    output_1300_1_perp.momZOut(cond_1300_1_perp,1),output_1300_1_perp.momZOut(cond_1300_1_perp,2),1,25);

KER_1300_par_1 = getCOMKER2(output_1300_1_par.momXOut(cond_1300_1_par,1),output_1300_1_par.momXOut(cond_1300_1_par,2),...
    output_1300_1_par.momYOut(cond_1300_1_par,1),output_1300_1_par.momYOut(cond_1300_1_par,2),...
    output_1300_1_par.momZOut(cond_1300_1_par,1),output_1300_1_par.momZOut(cond_1300_1_par,2),1,25);

KER_800_perp_1 = getCOMKER2(output_800_1_perp.momXOut(cond_800_1_perp,1),output_800_1_perp.momXOut(cond_800_1_perp,2),...
    output_800_1_perp.momYOut(cond_800_1_perp,1),output_800_1_perp.momYOut(cond_800_1_perp,2),...
    output_800_1_perp.momZOut(cond_800_1_perp,1),output_800_1_perp.momZOut(cond_800_1_perp,2),1,25);

KER_800_par_1 = getCOMKER2(output_800_1_par.momXOut(cond_800_1_par,1),output_800_1_par.momXOut(cond_800_1_par,2),...
    output_800_1_par.momYOut(cond_800_1_par,1),output_800_1_par.momYOut(cond_800_1_par,2),...
    output_800_1_par.momZOut(cond_800_1_par,1),output_800_1_par.momZOut(cond_800_1_par,2),1,25);

[hist266,X266] = hist(KER_266_1,histx);
[hist1300_perp,X1300] = hist(KER_1300_perp_1,histx);
[hist800_perp,X800] = hist(KER_800_perp_1,histx);
[hist1300_par,X1300] = hist(KER_1300_par_1,histx);
[hist800_par,X800] = hist(KER_800_par_1,histx);
figure()
errorbar(X266,hist266/2,sqrt(hist266)/2,sqrt(hist266)/2)
hold on;

norm_800_perp = hist800_perp/sum(hist800_perp);
error_800_perp = sqrt(hist800_perp)/sum(hist800_perp);
norm_800_par = hist800_par/sum(hist800_par);
error_800_par = sqrt(hist800_par)/sum(hist800_par);

errorbar(X800,hist800_perp,sqrt(hist800_perp),sqrt(hist800_perp))
errorbar(X800,hist800_par,sqrt(hist800_par),sqrt(hist800_par))
xlabel('KER (eV)')
legend('266','800 perp','800 par')
title('C2H/H 800')

error_266 = sqrt(hist266)/2;

figure()
errorbar(X266,hist266/2,sqrt(hist266)/2,sqrt(hist266)/2)
hold on;

norm_1300_perp = hist1300_perp/sum(hist1300_perp);
error_1300_perp = sqrt(hist1300_perp)/sum(hist1300_perp);
norm_1300_par = hist1300_par/sum(hist1300_par);
error_1300_par = sqrt(hist1300_par)/sum(hist1300_par);

errorbar(X1300,hist1300_perp,sqrt(hist1300_perp),sqrt(hist1300_perp))
errorbar(X1300,hist1300_par,sqrt(hist1300_par),sqrt(hist1300_par))
xlabel('KER (eV)')
legend('266','1300 perp','1300 par')
title('C2H/H 1300')
%%

figure()
errorbar(X800,(hist800_perp-hist266/2),sqrt(error_800_perp.^2+error_266.^2),sqrt(error_800_perp.^2+error_266.^2))
hold on
errorbar(X800,(hist800_par-hist266/2),sqrt(error_800_par.^2+error_266.^2),sqrt(error_800_par.^2+error_266.^2))
legend('800 perp', '800 par')
title('C2H/H 800')
xlabel('eV')
ylabel('change in counts from 266')


figure()
errorbar(X1300,(hist1300_perp-hist266/2),sqrt(error_1300_perp.^2+error_266.^2),sqrt(error_1300_perp.^2+error_266.^2))
hold on
errorbar(X1300,(hist1300_par-hist266/2),sqrt(error_1300_par.^2+error_266.^2),sqrt(error_1300_par.^2+error_266.^2))
legend('1300 perp', '1300 par')
title('C2H/H 1300')
xlabel('eV')
ylabel('change in counts from 266')

%%
BondMomX = output_266_1.momXOut(cond_266,1)-output_266_1.momXOut(cond_266,2);
BondMomY = output_266_1.momYOut(cond_266,1)-output_266_1.momYOut(cond_266,2);
BondMomZ = output_266_1.momZOut(cond_266,1)-output_266_1.momZOut(cond_266,2);

KER = output_266_1.KER(cond_266);

tt = '266';
sp = 'deprotonation';

KERAnglePlot


BondMomX = output_800_1_perp.momXOut(cond_800_1_perp,1)-output_800_1_perp.momXOut(cond_800_1_perp,2);
BondMomY = output_800_1_perp.momYOut(cond_800_1_perp,1)-output_800_1_perp.momYOut(cond_800_1_perp,2);
BondMomZ = output_800_1_perp.momZOut(cond_800_1_perp,1)-output_800_1_perp.momZOut(cond_800_1_perp,2);

KER = output_800_1_perp.KER(cond_800_1_perp);

tt = '800 perp';
sp = 'deprotonation';

KERAnglePlot


BondMomX = output_1300_1_perp.momXOut(cond_1300_1_perp,1)-output_1300_1_perp.momXOut(cond_1300_1_perp,2);
BondMomY = output_1300_1_perp.momYOut(cond_1300_1_perp,1)-output_1300_1_perp.momYOut(cond_1300_1_perp,2);
BondMomZ = output_1300_1_perp.momZOut(cond_1300_1_perp,1)-output_1300_1_perp.momZOut(cond_1300_1_perp,2);

KER = output_1300_1_perp.KER(cond_1300_1_perp);

tt = '1300 perp';
sp = 'deprotonation';

KERAnglePlot


BondMomX = output_800_1_par.momXOut(cond_800_1_par,1)-output_800_1_par.momXOut(cond_800_1_par,2);
BondMomY = output_800_1_par.momYOut(cond_800_1_par,1)-output_800_1_par.momYOut(cond_800_1_par,2);
BondMomZ = output_800_1_par.momZOut(cond_800_1_par,1)-output_800_1_par.momZOut(cond_800_1_par,2);

KER = output_800_1_par.KER(cond_800_1_par);

tt = '800 par';
sp = 'deprotonation';

KERAnglePlot


BondMomX = output_1300_1_par.momXOut(cond_1300_1_par,1)-output_1300_1_par.momXOut(cond_1300_1_par,2);
BondMomY = output_1300_1_par.momYOut(cond_1300_1_par,1)-output_1300_1_par.momYOut(cond_1300_1_par,2);
BondMomZ = output_1300_1_par.momZOut(cond_1300_1_par,1)-output_1300_1_par.momZOut(cond_1300_1_par,2);

KER = output_1300_1_par.KER(cond_1300_1_par);

tt = '1300 par';
sp = 'deprotonation';

KERAnglePlot
