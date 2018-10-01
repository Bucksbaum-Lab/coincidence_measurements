function Build_Fly2(charge, mass, posTOF, posX, posY, fly_spectromete_filename)
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

pos_strx = 'sequence {';
for ii = 1:length(posTOF)
    if ii == length(posTOF)
        pos_strx = horzcat(pos_strx,num2str(posTOF(ii)),'}');
    else
        pos_strx = horzcat(pos_strx,num2str(posTOF(ii)),', ');
    end
end

pos_stry = 'sequence {';
for ii = 1:length(posY)
    if ii == length(posY)
        pos_stry = horzcat(pos_stry,num2str(posY(ii)),'}');
    else
        pos_stry = horzcat(pos_stry,num2str(posY(ii)),', ');
    end
end

pos_strz = 'sequence {';
for ii = 1:length(posX)
    if ii == length(posX)
        pos_strz = horzcat(pos_strz,num2str(posX(ii)),'}');
    else
        pos_strz = horzcat(pos_strz,num2str(posX(ii)),', ');
    end
end

fly = 'particles {\n\tcoordinates = 0,\n';

fly = horzcat(fly,...
        '\tstandard_beam {\n',...
        '\t\ttob = 0,\n',...
        '\t\tmass = ',num2str(mass),',\n',...
        '\t\tcharge = 1,\n',...
        '\t\tcwf = 1,\n',...
        '\t\tcolor = 0,\n',...
        '\t\tx = ',num2str(pos_strx),',\n',...
        '\t\ty = ',num2str(pos_stry),',\n',...
        '\t\tz = ',num2str(pos_strz),',\n',...
        '\t\taz = 0,\n',...
        '\t\tel = 0,\n',...
        '\t\tke = 0\n');

fly = strcat(fly, '\t}\n}');

fileID = fopen([cd fly_spectromete_filename],'w+');

fprintf(fileID, fly);

fclose('all');
