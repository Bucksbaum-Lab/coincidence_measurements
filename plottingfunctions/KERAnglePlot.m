figure();

%BondPer = sqrt(BondMomX.^2+BondMomZ.^2);
%BondAngle = abs(atan(BondPer./BondMomY));
%BondAngle(BondAngle>pi()/2) = -BondAngle(BondAngle>pi()/2)+pi();
BondPer = sqrt(BondMomZ.^2+BondMomY.^2);
BondAngle = abs(atan(BondPer./BondMomX));

[N,C]=hist3([BondAngle,KER],[4,10]);

N=N';

Cx = C{1};
deltaX = Cx(2)-Cx(1);

Cy = C{2};
deltaY = Cy(2)-Cy(1);

X = [Cx,Cx(end)+deltaX];
Y = [Cy,Cy(end)+deltaY];
N = [N,zeros(length(Y)-1,1)];
N = [N;zeros(1,length(X))];

pcolor(X, Y, N);
title([tt ', ' sp])
xlabel('angle (rad)')
ylabel('KER (eV)')
colorbar
