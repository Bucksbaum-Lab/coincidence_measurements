%% Interpolate Over a Grid Using Cubic Method  

% Copyright 2015 The MathWorks, Inc.


%% 
% Coarsely sample the peaks function. 
[X,Y] = meshgrid(-3:3);
V = peaks(7);  

%% 
% Plot the coarse sampling. 
figure
surf(X,Y,V)
title('Original Sampling');     

%% 
% Create the query grid with spacing of 0.25. 
[Xq,Yq] = meshgrid(-3:0.25:3);  

%% 
% Interpolate at the query points, and specify cubic interpolation. 
Vq = interp2(X,Y,V,Xq,Yq,'cubic');  

%% 
% Plot the result. 
figure
surf(Xq,Yq,Vq);
title('Cubic Interpolation Over Finer Grid');      
