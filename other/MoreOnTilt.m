numx = 10;
numy = 10;

ay = 1;
ax = 2;
az = 1;

figure();
[outputy, centers] = hist3([ay*output_1300_1.momYOut(:,1), output_1300_1.partEnergyOut(:,1)], [numx,numy],'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputy)
title('Y')

figure();
[outputx, centers] = hist3([ax*output_1300_1.momXOut(:,1), output_1300_1.partEnergyOut(:,1)], [numx,numy],'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputx)
title('X')

figure();
[outputz, centers] = hist3([az*output_1300_1.momZOut(:,1), output_1300_1.partEnergyOut(:,1)], [numx,numy],'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputz)
title('Z')