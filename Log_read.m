function [Sim] = Log_read(file)
%Read SIMION log ASCII file for flight Sim and make into MATLAB matrix

    fid = fopen(file,'r+');
    fgetl(fid);
    buffer = fread(fid,Inf);
    fclose(fid);
    
    fid = fopen(file,'w+');
    fwrite(fid,buffer);
    
    fclose(fid);
    
    Sim = dlmread(file,',');

end

