%%
filename = 'G:\2018_05_11\take1\acetylene_1ps_6p4en10_torr_800_1300_266_DAn.txt';
nstrings = 10000000;
%%
fin  = fopen(filename);
fout = fopen([ filename(1:(length(filename)-4)), '_', num2str(nstrings, '%1.1e' ),...
               filename((length(filename)-3):length(filename)) ], 'w');

for i = 1:nstrings
    fprintf([fout], fgets(fin));
end

fclose(fin);
fclose(fout);