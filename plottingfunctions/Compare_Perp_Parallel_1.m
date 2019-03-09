histx = linspace(0,18,18*1+1);

C2HeV = 2;
HeV = 15;
mom = 7;
totKER = 18;

coinNum = linspace(1,size(output_266_1.numHitsOut,1),size(output_266_1.numHitsOut,1));
coinNum = coinNum';
cond_266 = abs(sum(output_266_1.momXOut,2))<mom&...
    abs(sum(output_266_1.momYOut,2))<mom&...
    abs(sum(output_266_1.momZOut,2))<mom&...
    output_266_1.KER<totKER&...
    output_266_1.partEnergyOut(:,1)<HeV&output_266_1.partEnergyOut(:,2)<C2HeV&...
... %    coinNum<max(coinNum)/2;
... %    coinNum>max(coinNum)/2;
    coinNum>=0;

coinNum = linspace(1,size(output_1300_1_perp.numHitsOut,1),size(output_1300_1_perp.numHitsOut,1));
coinNum = coinNum';
cond_1300_1_perp = abs(sum(output_1300_1_perp.momXOut,2))<mom&...
    abs(sum(output_1300_1_perp.momYOut,2))<mom&...
    abs(sum(output_1300_1_perp.momZOut,2))<mom&...
    output_1300_1_perp.KER<totKER&...
    output_1300_1_perp.partEnergyOut(:,1)<HeV&output_1300_1_perp.partEnergyOut(:,2)<C2HeV&...
... %    coinNum<max(coinNum)/2;
... %    coinNum>max(coinNum)/2;
    coinNum>=0;

coinNum = linspace(1,size(output_1300_1_par.numHitsOut,1),size(output_1300_1_par.numHitsOut,1));
coinNum = coinNum';
cond_1300_1_par = abs(sum(output_1300_1_par.momXOut,2))<mom&...
    abs(sum(output_1300_1_par.momYOut,2))<mom&...
    abs(sum(output_1300_1_par.momZOut,2))<mom&...
    output_1300_1_par.KER<totKER&...
    output_1300_1_par.partEnergyOut(:,1)<HeV&output_1300_1_par.partEnergyOut(:,2)<C2HeV&...
... %    coinNum<max(coinNum)/2;
... %    coinNum>max(coinNum)/2;
    coinNum>=0;

coinNum = linspace(1,size(output_800_1_perp.numHitsOut,1),size(output_800_1_perp.numHitsOut,1));
coinNum = coinNum';
cond_800_1_perp = abs(sum(output_800_1_perp.momXOut,2))<mom&...
    abs(sum(output_800_1_perp.momYOut,2))<mom&...
    abs(sum(output_800_1_perp.momZOut,2))<mom&...
    output_800_1_perp.KER<totKER&...
    output_800_1_perp.partEnergyOut(:,1)<HeV&output_800_1_perp.partEnergyOut(:,2)<C2HeV&...
... %    coinNum<max(coinNum)/2;
... %    coinNum>max(coinNum)/2;
    coinNum>=0;

coinNum = linspace(1,size(output_800_1_par.numHitsOut,1),size(output_800_1_par.numHitsOut,1));
coinNum = coinNum';
cond_800_1_par = abs(sum(output_800_1_par.momXOut,2))<mom&...
    abs(sum(output_800_1_par.momYOut,2))<mom&...
    abs(sum(output_800_1_par.momZOut,2))<mom&...
    output_800_1_par.KER<totKER&...
    output_800_1_par.partEnergyOut(:,1)<HeV&output_800_1_par.partEnergyOut(:,2)<C2HeV&...
... %    coinNum<max(coinNum)/2;
... %    coinNum>max(coinNum)/2;
    coinNum>=0;

KER_266_1 = output_266_1.KER(cond_266);
KER_1300_perp_1 = output_1300_1_perp.KER(cond_1300_1_perp);
KER_800_perp_1 = output_800_1_perp.KER(cond_800_1_perp);
KER_800_par_1 = output_800_1_par.KER(cond_800_1_par);
KER_1300_par_1 = output_1300_1_par.KER(cond_1300_1_par);

[hist266,X266] = hist(KER_266_1,histx);
[hist1300_perp,X1300] = hist(KER_1300_perp_1,histx);
[hist800_perp,X800] = hist(KER_800_perp_1,histx);
[hist1300_par,X1300] = hist(KER_1300_par_1,histx);
[hist800_par,X800] = hist(KER_800_par_1,histx);

error_266 = sqrt(hist266)/2;

error_800_perp = sqrt(hist800_perp);
error_800_par = sqrt(hist800_par);

error_1300_perp = sqrt(hist1300_perp);
error_1300_par = sqrt(hist1300_par);
%{
figure('units', 'inches', 'position', [.5 .5 7.5 6.5])
errorbar(X266,hist266/2,error_266,error_266,'linewidth', 3, 'color',[0.4660, 0.6740, 0.1880])
hold on;

errorbar(X800,hist800_perp,error_800_perp,error_800_perp,'linewidth', 3, 'color', [0,112,184]/255)
errorbar(X800,hist800_par,error_800_par,error_800_par,'linewidth', 3, 'color', 'k')
xlabel('KER (eV)')
ylabel('counts')
legend('266','800 perp','800 par')
title('Deprotonation, 800 nm Control')
set(gca, 'fontsize', 18)
xlim([0,18])

figure('units', 'inches', 'position', [.5 .5 7.5 6.5])
errorbar(X266,hist266/2,error_266,error_266,'linewidth', 3, 'color',[0.4660, 0.6740, 0.1880])
hold on;

errorbar(X1300,hist1300_perp,error_1300_perp,error_1300_perp,'linewidth', 3, 'color', [0,112,184]/255)
errorbar(X1300,hist1300_par,error_1300_par,error_1300_par,'linewidth', 3, 'color', 'k')
xlabel('KER (eV)')
ylabel('counts')
legend('266','1300 perp','1300 par')
title('Deprotonation, 1300 nm Control')
set(gca, 'fontsize', 18)
xlim([0,18])


figure('units', 'inches', 'position', [.5 .5 7.5 6.5])
errorbar(X800,(hist800_perp-hist266/2),sqrt(error_800_perp.^2+error_266.^2),sqrt(error_800_perp.^2+error_266.^2),'linewidth', 3', 'color', [0,112,184]/255)
hold on
errorbar(X800,(hist800_par-hist266/2),sqrt(error_800_par.^2+error_266.^2),sqrt(error_800_par.^2+error_266.^2),'linewidth', 3, 'color', 'k')
xlim([0,18])
legend('perp', 'par')
title('Deprotonation, 800 nm Control')
xlabel('KER (eV)')
ylabel('\Delta counts')
set(gca, 'fontsize', 18)


figure('units', 'inches', 'position', [.5 .5 7.5 6.5])
errorbar(X1300,(hist1300_perp-hist266/2),sqrt(error_1300_perp.^2+error_266.^2),sqrt(error_1300_perp.^2+error_266.^2),'linewidth', 3', 'color', [0,112,184]/255)
hold on
errorbar(X1300,(hist1300_par-hist266/2),sqrt(error_1300_par.^2+error_266.^2),sqrt(error_1300_par.^2+error_266.^2),'linewidth', 3, 'color', 'k')
xlim([0,18])
legend('perp', 'par')
title('Deprotonation, 1300 nm Control')
xlabel('KER (eV)')
ylabel('\Delta counts')
set(gca, 'fontsize', 18)
%}

figure('units', 'inches', 'position', [.5 .5 7.5 5.5])
subplot(2,1,1)
errorbar(X1300,hist1300_par,error_1300_par,error_1300_par,'-*','linewidth', 2, 'color', [0,112,184]/255)
hold on
errorbar(X1300,hist1300_perp,error_1300_perp,error_1300_perp,'linewidth', 2, 'color', 'k')
errorbar(X266,hist266/2,error_266,error_266,'--*','linewidth', 2, 'color', [0.4660, 0.6740, 0.1880])
hold off
xlim([0,18])
legend('par','perp', 'off')
title('1300 nm control')
xlabel('KER (eV)')
ylabel('counts')
subplot(2,1,2)
errorbar(X800,hist800_perp,error_800_perp,error_800_perp,'-*','linewidth', 2, 'color', [0,112,184]/255)
hold on
errorbar(X800,hist800_par,error_800_par,error_800_par,'linewidth', 2, 'color', 'k')
errorbar(X266,hist266/2,error_266,error_266,'--*','linewidth', 2, 'color', [0.4660, 0.6740, 0.1880])
xlim([0,18])
legend('par', 'perp', 'off')
title('800 nm control')
xlabel('KER (eV)')
ylabel('counts')

%[.5 .5 4.5 3.5]
figure('units', 'inches', 'position', [.5 .5 7.5 5.5])
subplot(2,1,1)
errorbar(X1300,(hist1300_par-hist266/2),sqrt(error_1300_par.^2+error_266.^2),sqrt(error_1300_par.^2+error_266.^2),'-*','linewidth', 2, 'color', [0,112,184]/255)
hold on
errorbar(X1300,(hist1300_perp-hist266/2),sqrt(error_1300_perp.^2+error_266.^2),sqrt(error_1300_perp.^2+error_266.^2),'linewidth', 2, 'color', 'k')
hold off
xlim([0,18])
legend('par', 'perp')
title('1300 nm control')
xlabel('KER (eV)')
ylabel('\Delta counts')
subplot(2,1,2)
errorbar(X800,(hist800_par-hist266/2),sqrt(error_800_par.^2+error_266.^2),sqrt(error_800_par.^2+error_266.^2),'-*','linewidth', 2, 'color', [0,112,184]/255)
hold on
errorbar(X800,(hist800_perp-hist266/2),sqrt(error_800_perp.^2+error_266.^2),sqrt(error_800_perp.^2+error_266.^2),'linewidth', 2, 'color', 'k')
xlim([0,18])
legend('par', 'perp')
title('800 nm control')
xlabel('KER (eV)')
ylabel('\Delta counts')
