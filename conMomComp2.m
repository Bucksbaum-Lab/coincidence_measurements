function cond = conMomComp2(mv1z, mv2z, mv1x, mv2x, mv1y, mv2y, deltat, deltax, deltay)

if deltat > 0
    cond = mv1z+mv2z < deltat & mv1z+mv2z > -deltat;
else
    cond = true(size(mv1x));
end

if deltax > 0
    cond = cond & mv1x+mv2x < deltax & mv1x+mv2x > -deltax;
end

if deltay > 0
    cond = cond & mv1y+mv2y < deltay & mv1y+mv2y > -deltay;
end
