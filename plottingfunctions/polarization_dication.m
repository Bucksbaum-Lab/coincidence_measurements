close all

figure(); hist(output_266_13.momXOut(:),20)
title('13,X')
figure(); hist(output_266_13.momYOut(:),20)
title('13,Y')
figure(); hist(output_266_13.momZOut(:),20)
title('13,Z')
figure(); hist(output_266_14.momXOut(:),20)
title('14,X')
figure(); hist(output_266_14.momYOut(:),20)
title('14,Y')
figure(); hist(output_266_14.momZOut(:),20)
title('14,Z')

Centers = {-15:2:15 0:2:15};

KER = repmat(output_266_14.KER,1,2);
figure();
hist3([output_266_14.momYOut(:),KER(:)],'CdataMode', 'auto','Ctrs',Centers)
xlabel('Ymom')
ylabel('KER')
title('14')
view(2)
colorbar

KER = repmat(output_266_13.KER,1,2);
figure();
hist3([output_266_13.momYOut(:),KER(:)],'CdataMode', 'auto','Ctrs',Centers)
xlabel('Ymom')
ylabel('KER')
title('13')
view(2)
colorbar

KER = repmat(output_266_14.KER,1,2);
figure();
hist3([output_266_14.momXOut(:),KER(:)],'CdataMode', 'auto','Ctrs',Centers)
xlabel('Xmom')
ylabel('KER')
title('14')
view(2)
colorbar

KER = repmat(output_266_13.KER,1,2);
figure();
hist3([output_266_13.momXOut(:),KER(:)],'CdataMode', 'auto','Ctrs',Centers)
xlabel('Xmom')
ylabel('KER')
title('13')
view(2)
colorbar

KER = repmat(output_266_14.KER,1,2);
figure();
hist3([output_266_14.momZOut(:),KER(:)],'CdataMode', 'auto','Ctrs',Centers)
xlabel('Zmom')
ylabel('KER')
title('14')
view(2)
colorbar

KER = repmat(output_266_13.KER,1,2);
figure();
hist3([output_266_13.momZOut(:),KER(:)],'CdataMode', 'auto','Ctrs',Centers)
xlabel('Zmom')
ylabel('KER')
title('13')
view(2)
colorbar










%%

CCaxis_13 = [output_266_13.momXOut(:,1)-output_266_13.momXOut(:,2),...
    output_266_13.momYOut(:,1)-output_266_13.momYOut(:,2),...
    output_266_13.momZOut(:,1)-output_266_13.momZOut(:,2)];
CCaxis_14 = [output_266_14.momXOut(:,1)-output_266_14.momXOut(:,2),...
    output_266_14.momYOut(:,1)-output_266_14.momYOut(:,2),...
    output_266_14.momZOut(:,1)-output_266_14.momZOut(:,2)];

CosAngle_13 = acos(abs(CCaxis_13(:,2)./sqrt((CCaxis_13(:,1)).^2+(CCaxis_13(:,2)).^2+(CCaxis_13(:,3)).^2)));

CosAngle_14 = acos(abs(CCaxis_14(:,2)./sqrt((CCaxis_14(:,1)).^2+(CCaxis_14(:,2)).^2+(CCaxis_14(:,3)).^2)));

Centers = {0:pi/5:pi/2 0:1:10};
figure();
hist3([abs(CosAngle_13),output_266_13.KER],'CdataMode', 'auto', 'Ctrs', Centers);
title('13')
colorbar
view(2)

figure();
hist3([abs(CosAngle_14),output_266_14.KER],'CdataMode', 'auto', 'Ctrs', Centers);
title('14')
colorbar
view(2)

hist13 = hist(abs(CosAngle_13),0:pi/20:pi/2);
hist14 = hist(abs(CosAngle_14),0:pi/20:pi/2);

figure();
plot(0:pi/20:pi/2, hist13/sum(hist13), '-o')
hold on
plot(0:pi/20:pi/2, hist14/sum(hist14), '-o')
legend('13', '14')
















