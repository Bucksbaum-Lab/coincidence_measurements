x = 0;

for ii = 1:192
    x=[x;1;1;0;0;1;1;0;1;0;1;0;1;0];
end

save('ShutterClosed.txt', 'x', '-ascii')