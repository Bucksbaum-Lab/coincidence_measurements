function Locations = get2peaks(Y)

[~,loc] = findpeaks(Y,'SortStr','descend');
try
    Locations = [loc(1), loc(2)];
catch
    Locations = [loc(1), loc(1)];
end

end