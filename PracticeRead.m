ds = tabularTextDatastore('10000points.txt');
numchunks = 3;
EventTagn = [EventTag, 10^9];
for nn=1:(length(EventTagn)-1)/3
    
    ds.ReadSize = EventTagn(numchunks*nn+1)-EventTagn(numchunks*nn-(numchunks-1));
    T = read(ds);
    
end
