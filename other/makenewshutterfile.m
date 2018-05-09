x = 0;

for ii = 1:164
    x=[x;1;1;0;0;1;1;0;1;0;1;0;1;0];
end

x = [x;1];
save('ShutterClosed.txt', 'x', '-ascii')