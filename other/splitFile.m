%%
filename = 'G:\2018_03_09\Acetylene_5p4en10torr_266_5p6mW_1300_38mW_800_103mW_500fsDelay_DAn_DAn.txt';
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