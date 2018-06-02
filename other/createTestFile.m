initialPosition = 37;
V1 = 2000;
VM = -2250;
testMass1 = 1;
testCharge1 = 1;
testMass2 = 12;
testCharge2 = 1;
energies = linspace(0,10,11);
thetas = linspace(0,360,10);
phis = linspace(0,360,10);

%retreive fly'm data from simion
Sim = Flym_Sim(testCharge1, testMass1, energies, thetas, phis, initialPosition, V1, VM);
t_Sim = Sim(2:2:end, 2)*10^3;
x_Sim = Sim(2:2:end, 3)-55;
y_Sim = Sim(2:2:end, 4)-55;

Output = zeros(length(t_Sim),14);

Output(:,2) = 3;
Output(:,1) = linspace(1,length(t_Sim),length(t_Sim));
Output(:,3) = x_Sim;
Output(:,4) = y_Sim;
Output(:,5) = t_Sim;


%retreive fly'm data from simion
Sim = Flym_Sim(testCharge1, testMass1, energies, thetas-180, -phis, initialPosition, V1, VM);
t_Sim = Sim(2:2:end, 2)*10^3;
x_Sim = Sim(2:2:end, 3)-55;
y_Sim = Sim(2:2:end, 4)-55;

Output(:,6) = x_Sim;
Output(:,7) = y_Sim;
Output(:,8) = t_Sim;
%}
%retreive fly'm data from simion
Sim = Flym_Sim(testCharge2, testMass2, energies, thetas, phis, initialPosition, V1, VM);
t_Sim = Sim(2:2:end, 2)*10^3;
x_Sim = Sim(2:2:end, 3)-55;
y_Sim = Sim(2:2:end, 4)-55;

Output(:,9) = x_Sim;
Output(:,10) = y_Sim;
Output(:,11) = t_Sim;

dlmwrite('fakeData.txt', Output, 'delimiter', ' ', 'newline', 'pc')

