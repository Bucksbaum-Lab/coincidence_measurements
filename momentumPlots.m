clear all


dpx = 3;
filename= 'G:\\2018_02_09\\analysis\\Acetylene_266_800_1300_4p3e9torr_80fs_DAn_DAn-cut-23-Feb-2018-22-4-41_M=1,12,13_Z=1(10eV, dpx=3)_OPEN.mat';
load(filename);

figure()
subplot(1,3,1) 
hist(sum(real(output.momXOut),2), 50)
subplot(1,3,2)
hist(sum(real(output.momYOut),2), 50)
subplot(1,3,3)
hist(sum(real(output.momZOut),2), 50)

momX = output.momXOut;
momY = output.momYOut;
momZ = output.momZOut;

mass = output.mass;

p = 1;
p_on = [2,3];

KER = output.KER;
    
CMCX = (momX(:, p_on(1))+momX(:, p_on(2)))/sum(mass(p_on));
CMCY = (momY(:, p_on(1))+momY(:, p_on(2)))/sum(mass(p_on));
CMCZ = (momZ(:, p_on(1))+momZ(:, p_on(2)))/sum(mass(p_on));


momX(:, p) = momX(:, p) - CMCX*mass(p);
momX(:, p_on(1)) = momX(:, p_on(1)) - CMCX*mass(p_on(1));
momX(:, p_on(2)) = momX(:, p_on(2)) - CMCX*mass(p_on(2));

momY(:, p) = momY(:, p) - CMCX*mass(p);
momY(:, p_on(1)) = momY(:, p_on(1)) - CMCX*mass(p_on(1));
momY(:, p_on(2)) = momY(:, p_on(2)) - CMCX*mass(p_on(2));

momZ(:, p) = momZ(:, p) - CMCX*mass(p);
momZ(:, p_on(1)) = momZ(:, p_on(1)) - CMCX*mass(p_on(1));
momZ(:, p_on(2)) = momZ(:, p_on(2)) - CMCX*mass(p_on(2));

V = [(momX(:, p_on(2))-momX(:, p_on(1))), ...
     (momY(:, p_on(2))-momY(:, p_on(1))), ...
     (momZ(:, p_on(2))-momZ(:, p_on(1)))];
V = V./apply_to_rows(@norm, V);

parallel_proj  = dot(V, [momX(:, p), momY(:, p), momZ(:, p)], 2);
perpendicular_proj = apply_to_rows(@norm, [momX(:, p), momY(:, p), momZ(:, p)] - V.*parallel_proj); 
%%
figure()
[N,C] = hist3([(abs(real(parallel_proj))), real(perpendicular_proj)], [30, 30]);
pcolor(C{1}, C{2}, N');
title('shutter open')
xlabel('$$ \vec{P}(H^{+})\parallel \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')
ylabel('$$ \vec{P}(H^{+})\perp \left(\vec{P}(CH^{+}) - \vec{P}(C^{+})\right) $$', 'Interpreter','latex')