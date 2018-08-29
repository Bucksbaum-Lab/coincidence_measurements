points = load('C:\Data\testdata\100000points.txt');

pointsofinterest = points(:,5)>2640 & points(:,5)<2670;
xpoints = points(pointsofinterest,3);
ypoints = points(pointsofinterest,4);
inverset = points(pointsofinterest,5);
inverset = 1./inverset;

xpoints = xpoints-median(xpoints);
ypoints = ypoints-median(ypoints);

pointsofinterest = (xpoints < 2.5)&(xpoints > -2.5)&(ypoints < 2.5)&(ypoints > -2.5);

xpoints = xpoints(pointsofinterest);
ypoints = ypoints(pointsofinterest);
inverset = inverset(pointsofinterest);

xpoints = xpoints-mean(xpoints);
ypoints = ypoints-mean(ypoints);
inverset = inverset - mean(inverset);

xpoints = xpoints/std(xpoints);
ypoints = ypoints/std(ypoints);
inverset = inverset/std(inverset);


numx = 20;
numy = 20;

figure();
[outputy, centers] = hist3([ypoints,inverset], [numx,numy],'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputy)

hold on
[outputx, centers] = hist3([xpoints,inverset], centers,'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputx,'--')
title('tof vs x or y')
legend('tof v y', 'tof v x')
view(2)
ylim([-2.5 2.5])
xlim([-2.5 2.5])

figure();
contour(X',Y',outputy-outputx)
title('tof vs y - tof vs x')
view(2)
ylim([-2.5 2.5])
xlim([-2.5 2.5])






figure();
[outputy, centers] = hist3([inverset, xpoints], [numx,numy],'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputy)

hold on
[outputx, centers] = hist3([xpoints,inverset], [numx, numy],'CdataMode', 'auto');
x = centers{:,1};
y = centers{:,2};
[X,Y] = meshgrid(x,y);
contour(X',Y',outputx,'--')
title('x vs tof or y')
legend('x v tof', 'x v tof')
view(2)
ylim([-2.5 2.5])
xlim([-2.5 2.5])

figure();
contour(X',Y',outputy-outputx)
title('x vs tof - x vs tof')
view(2)
ylim([-2.5 2.5])
xlim([-2.5 2.5])