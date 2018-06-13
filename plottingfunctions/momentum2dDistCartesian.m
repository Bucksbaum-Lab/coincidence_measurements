function[X, Y, N] = momentum2dDistCartesian ...
            (parallel_proj, perpendicular_proj, nbins, plotname, name1, name2, name3)

    [N,C] = hist3([real(parallel_proj), real(perpendicular_proj)], [nbins, nbins]);
    N = N';
    
    Xp = C{1}(1:(length(C{1}) - 1)) + diff(C{1})/2 ;
    Xm = C{1}(1:(length(C{1}) - 1)) - diff(C{1})/2 ;
    
    Yp = C{2}(1:(length(C{2}) - 1)) + diff(C{2})/2 ;
    Ym = C{2}(1:(length(C{2}) - 1)) - diff(C{2})/2 ;
    
    X = [Xm(1), (Xm(2:length(Xm)) + Xp(1:(length(Xp) - 1)))/2, Xp(length(Xp))];
    Y = [Ym(1), (Ym(2:length(Xm)) + Yp(1:(length(Xp) - 1)))/2, Yp(length(Yp))];
    pcolor(X, Y, N);
    
    title(plotname)
    xlabel(['$$ \vec{P}(' name1 ')\parallel \left(\vec{P}(' name3 ') - \vec{P}(' name2 ')\right) $$'], 'Interpreter','latex')
    ylabel(['$$ \vec{P}(' name1 ')\perp  \left(\vec{P}('    name3 ') - \vec{P}(' name2 ')\right) $$'], 'Interpreter','latex')