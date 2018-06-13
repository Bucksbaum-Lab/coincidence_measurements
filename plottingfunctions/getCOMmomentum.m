function output = getCOMmomentum(m1,p1,m2,p2,m3,p3)

V = (p1+p2+p3)/(m1+m2+m3);

P1 = p1-m1*V;
P2 = p2-m2*V;
P3 = p3-m3*V;

output = [P1, P2, P3];