% Xin_full =  normrnd(2, 3, 100000000, 1);
% Yin_full =  normrnd(0, 10, 100000000, 1);   
% 
% Xin = Xin_full(Xin_full.^2 + Yin_full.^2 < 35^2);
% Yin = Yin_full(Xin_full.^2 + Yin_full.^2 < 35^2);
% 
% step = pi/64;
% rdata(2) = 0.;
% rdata(3) = 35;
% rdata(1) = 512;
% 
% %%
% H = hist3polar(Xin, Yin, 512, 0, 35, 128, -pi/2, pi/2);
% %%
% 
% pl = pcolor(H{1}, H{2}, H{3});
% set(pl, 'EdgeColor', 'none');
% pl
% axis equal tight


function [X, Y, Nout, rs] = hist3polar (Xin, Yin, rdata, adata)

    r_in = sqrt(Xin.^2 + Yin.^2);

    theta_in = atan(Yin./Xin) + (Xin < 0 & Yin >= 0)*pi + 2*(Xin < 0 & Yin < 0)*pi;
   
    rs = sqrt(rdata(2)^2:((rdata(3)^2 - rdata(2)^2))/(rdata(1)):rdata(3)^2);

    angles = adata(2):((adata(3)-adata(2))/adata(1)):(adata(3));

    [Nout, Cout] = hist3([r_in, theta_in], 'Edges', {rs, angles} );

    X = rs'.*cos(angles);
    Y = rs'.*sin(angles);
 end

