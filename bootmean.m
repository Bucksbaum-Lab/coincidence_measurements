figure()
[X,Y] = hist(bootstrp(1000,@(x) meanmean(x),plotting_data(1).polarDistNorm),10);
plot(Y,X,'-o')
hold on
[X,Y] = hist(bootstrp(1000,@(x) meanmean(x),plotting_data(2).polarDistNorm),10);
plot(Y,X,'-o')
[X,Y] = hist(bootstrp(1000,@(x) meanmean(x),plotting_data(3).polarDistNorm),10);
plot(Y,X,'-o')
[X,Y] = hist(bootstrp(1000,@(x) meanmean(x),plotting_data(4).polarDistNorm),10);
plot(Y,X,'-o')
legend('1','2','3','4')