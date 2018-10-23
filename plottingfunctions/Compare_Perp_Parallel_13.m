histx = linspace(0,18,18*2+1);

CHeV = 15;
mom = 10;
totKER = 18;

cond_266 = sum(output_266_13.momXOut,2)<mom&output_266_13.KER<totKER&...
    output_266_13.partEnergyOut(:,1)<CHeV&output_266_13.partEnergyOut(:,2)<CHeV;
cond_1300_13_perp = sum(output_1300_13_perp.momXOut,2)<mom&output_1300_13_perp.KER<totKER&...
    output_1300_13_perp.partEnergyOut(:,1)<CHeV&output_1300_13_perp.partEnergyOut(:,2)<CHeV;
cond_1300_13_par = sum(output_1300_13_par.momXOut,2)<mom&output_1300_13_par.KER<totKER&...
    output_1300_13_par.partEnergyOut(:,1)<CHeV&output_1300_13_par.partEnergyOut(:,2)<CHeV;
cond_800_13_perp = sum(output_800_13_perp.momXOut,2)<mom&output_800_13_perp.KER<totKER&...
    output_800_13_perp.partEnergyOut(:,1)<CHeV&output_800_13_perp.partEnergyOut(:,2)<CHeV;
cond_800_13_par = sum(output_800_13_par.momXOut,2)<mom&output_800_13_par.KER<totKER&...
    output_800_13_par.partEnergyOut(:,1)<CHeV&output_800_13_par.partEnergyOut(:,2)<CHeV;

[hist266,X266] = hist(output_266_13.KER(cond_266),histx);
[hist1300_perp,X1300] = hist(output_1300_13_perp.KER(cond_1300_13_perp),histx);
[hist800_perp,X800] = hist(output_800_13_perp.KER(cond_800_13_perp),histx);
[hist1300_par,X1300] = hist(output_1300_13_par.KER(cond_1300_13_par),histx);
[hist800_par,X800] = hist(output_800_13_par.KER(cond_800_13_par),histx);

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

figure()
errorbar(X800,(hist800_perp-hist800_par),sqrt(error_800_perp.^2+error_800_par.^2),sqrt(error_800_perp.^2+error_800_par.^2))
hold on
errorbar(X1300,(hist1300_perp-hist1300_par),sqrt(error_1300_perp.^2+error_1300_par.^2),sqrt(error_1300_perp.^2+error_1300_par.^2))
legend('800', '1300')
title('perp-par')
xlabel('eV')
ylabel('change in counts')