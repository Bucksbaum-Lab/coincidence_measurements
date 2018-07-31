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