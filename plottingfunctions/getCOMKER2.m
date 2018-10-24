function KERc = getCOMKER2(p1x,p2x,p1y,p2y,p1z,p2z,m1,m2)

Vcx = (p1x+p2x)/(m1+m2);

p1xc = p1x-m1*Vcx;
p2xc = p2x-m2*Vcx;


Vcy = (p1y+p2y)/(m1+m2);

p1yc = p1y-m1*Vcy;
p2yc = p2y-m2*Vcy;


Vcz = (p1z+p2z)/(m1+m2);

p1zc = p1z-m1*Vcz;
p2zc = p2z-m2*Vcz;



KERc = .5*(p1xc.^2/m1+p2xc.^2/m2 + p1yc.^2/m1+p2yc.^2/m2 + p1zc.^2/m1+p2zc.^2/m2);