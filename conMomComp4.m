function [cond] = conMomComp4(mv1, mv2, mv3, mv4, deltamv)

if deltamv > 0
    cond = mv1+mv2+mv3+mv4 < deltamv & mv1+mv2+mv3+mv4 > -deltamv;
else
    cond = true(size(mv1));
end
