function [Sim] = Flym_Sim(charge, mass, eV, theta, phi, ss, V1, VM, uniqueid)

%phi(eV < 0) = phi(eV < 0) + 180;
fly_spectromete_filename = ['.\flym\fly_Spectrometer'...
                            '_uniqueid=' uniqueid '.txt' '.fly2'];
record_out_filename = ['.\flym\Lab_record_out'...
                       '_uniqueid=' uniqueid '.txt'];

Build_Fly(charge, mass, eV, theta, phi, 53-ss, fly_spectromete_filename);

%delete file old record file
if exist(record_out_filename, 'file')
    delete(record_out_filename); 
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
Sim = Log_read(record_out_filename);%([cd '\flym\Lab_record_out.txt']);

end