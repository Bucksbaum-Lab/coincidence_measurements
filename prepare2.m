%{
%--------------------------------prepare----------------------------------%
%-------------------------------Version 4.0-------------------------------%

prepare generates matrices that will be used by cutandplot.  It generates
these matricies from position and tof data from LCLS data.  It will
generate matrices for X, Y, and Z with column representing the momentum
assuming a certain input mass.  It will also generate a vector that
indicates which hit and shot each entry comes from.

All of this is nicely opperated via the 'main.m' file which is a GUI.

%---------------------Last Edited on November 1, 2017---------------------%
%------Last Edited by and Original Author - Chelsea Liekhus-Schmaltz------% 
%}

function [EV, mom_tof, mom_x, mom_y] =...
    prepare2(V1, VM, ss, t0, x0, y0,...
             mass, charge, maxEV, ...
             tof, rX, rY, numHits, EVlength, Thetalength)

    [numshots, maxions] = size(tof);

    %initialize the momentum matrices
    EV      = cell(maxions,1);
    mom_tof = cell(maxions,1);
    mom_x   = cell(maxions,1);
    mom_y   = cell(maxions,1);

    for nioin = 1:maxions
        EV{nioin}      = sparse(length(rX), length(mass));
        mom_tof{nioin} = sparse(length(rX), length(mass));
        mom_x{nioin}   = sparse(length(rX), length(mass));
        mom_y{nioin}   = sparse(length(rX), length(mass));
    end
    %find the momentum associated with each event assuming it has a given mass
    %determined by which column it is in
    parfor n = 1:maxions
        cond = (numHits >= n);
        %parfor
        for ii = 1:length(mass)
            [EV{n}(cond, ii), mom_tof{n}(cond, ii), mom_x{n}(cond, ii), mom_y{n}(cond, ii)] =...
                convertToEnergy(tof(cond, n) - t0, rX(cond, n) - x0, rY(cond, n) - y0, ...
                                V1, VM, ss, charge(ii), mass(ii), maxEV(ii),...
                                EVlength, Thetalength, num2str(n));
        end
    end
end