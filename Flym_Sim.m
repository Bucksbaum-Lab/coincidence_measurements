function [Sim] = Flym_Sim(charge, mass, eV, theta, phi, ss, V1, VM)

%phi(eV < 0) = phi(eV < 0) + 180;

Build_Fly(charge, mass, eV, theta, phi, 53-ss);

%delete file old record file
if exist([cd '\flym\Lab_record_out.txt'], 'file')
    delete([cd '\flym\Lab_record_out.txt']); 
end

%run fly file in simion
system(['simion --nogui fastadj ' cd '\flym\TopGem.pa# 1=' num2str(V1)...
    ',2=' num2str(V1*2/3) ',3=' num2str(V1/3) ',4=0,5=' num2str(VM) ',6=0']);
system(['simion --nogui fly --retain-trajectories=0 --restore-potentials=0 --particles=' cd...
    '\flym\fly_Spectrometer.fly2 --recording-output=' cd '\flym\Lab_record_out.txt --recording='...
    cd '\flym\Lab_Spectrometer.rec ' cd '\flym\Lab_Spectrometer.iob']);

currentPath = cd;
cd flym
delete *.tmp
cd(currentPath)

%retreive flym data and display tof
Sim = Log_read([cd '\flym\Lab_record_out.txt']);

end