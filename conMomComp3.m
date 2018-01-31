function cond = conMomComp2(mv1, mv2, mv3, deltamv)

if deltamv > 0
    cond = (mv1+mv2+mv3) < deltamv & (mv1+mv2+mv3) > -deltamv;
else
    cond = true(size(mv1));
end
