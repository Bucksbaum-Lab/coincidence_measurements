function Build_Fly(charge, mass, energy, theta, phi, posTOF, fly_spectromete_filename)
    %Make fly2 file for determing tk for given k and n for different
        %transverse momentum

    %E is array of KEs. Each of these will be projected at various angles
        %which distributes the momentum between z, x, and y
    %Theta is the array of elevation angles (between pz and py) at which the Es are distributed over.
        %Use between 0 and 90 degrees 
        %This is elevation angle in SIMION
    %Phi is array of azimuthal angles (between px and pz)
        %0 is forward along TOF axis, 180 is backwards along TOF axis

mass = mass./charge;

E_str = 'sequence {';
for ii = 1:length(energy)
    if ii == length(energy)
        E_str = horzcat(E_str,num2str(energy(ii)),'}');
    else
        E_str = horzcat(E_str,num2str(energy(ii)),',');
    end
end

phi_str = 'sequence {';
for ii = 1:length(phi)
    if ii == length(phi)
        phi_str = horzcat(phi_str,num2str(phi(ii)),'}');
    else
        phi_str = horzcat(phi_str,num2str(phi(ii)),',');
    end
end

theta_str = 'sequence {';
for ii = 1:length(theta)
    if ii == length(theta)
        theta_str = horzcat(theta_str,num2str(theta(ii)),'}');
    else
        theta_str = horzcat(theta_str,num2str(theta(ii)),',');
    end
end

mass_str = 'sequence {';
for ii = 1:length(mass)
    if ii == length(mass)
        mass_str = horzcat(mass_str,num2str(mass(ii)),'}');
    else
        mass_str = horzcat(mass_str,num2str(mass(ii)),',');
    end
end

pos_str = 'sequence {';
for ii = 1:length(posTOF)
    if ii == length(posTOF)
        pos_str = horzcat(pos_str,num2str(posTOF(ii)),'}');
    else
        pos_str = horzcat(pos_str,num2str(posTOF(ii)),', ');
    end
end

fly = 'particles {\n\tcoordinates = 0,\n';

fly = horzcat(fly,...
        '\tstandard_beam {\n',...
        '\t\ttob = 0,\n',...
        '\t\tmass = ',mass_str,',\n',...
        '\t\tcharge = 1,\n',...
        '\t\tcwf = 1,\n',...
        '\t\tcolor = 0,\n',...
        '\t\tx = ',num2str(pos_str),',\n',...
        '\t\ty = 55,\n',...
        '\t\tz = 55,\n',...
        '\t\taz = ',theta_str,',\n',...
        '\t\tel = ',phi_str,',\n',...
        '\t\tke = ',E_str,'\n');

fly = strcat(fly, '\t}\n}');

fileID = fopen([cd fly_spectromete_filename],'w+');

fprintf(fileID, fly);

fclose('all');
