function [Sim] = Flym_Sim2(charge, mass, ss, x0, y0, V1, VM)

%phi(eV < 0) = phi(eV < 0) + 180;
fly_spectromete_filename = ['\flym\fly_Spectrometer'...
                            '_M=' strrep(num2str(mass),   '  ', '_') ...
                            '_Z=' strrep(num2str(charge), '  ', '_') '.fly2'];
record_out_filename = ['\flym\Lab_record_out'...
                       '_M=' strrep(num2str(mass),   '  ', '_') ...
                       '_Z=' strrep(num2str(charge), '  ', '_') '.txt'];

if exist(['.' fly_spectromete_filename], 'file')
    delete(['.' fly_spectromete_filename]);
end
                   
Build_Fly2(charge, mass, ss, x0, y0, fly_spectromete_filename);

%delete file old record file
if exist(['.' record_out_filename], 'file')
    delete(['.' record_out_filename]); 
end

%run fly file in simion
system(['simion --nogui fastadj ' cd '\flym\TopGem.pa# 1=' num2str(V1)...
    ',2=' num2str(V1*2/3) ',3=' num2str(V1/3) ',4=0,5=' num2str(VM) ',6=0']);

system(['simion --nogui fly --retain-trajectories=0 --restore-potentials=0' ...
    [' --particles=' [cd fly_spectromete_filename]] ...
    [' --recording-output=' [cd record_out_filename]] ...
    [' --recording=' [cd '\flym\Lab_Spectrometer.rec'] ' ' [cd '\flym\Lab_Spectrometer.iob']] ]);

currentPath = cd;
cd flym
delete *.tmp
cd(currentPath)

%retreive flym data and display tof
Sim = Log_read2(record_out_filename);%([cd '\flym\Lab_record_out.txt']);

end