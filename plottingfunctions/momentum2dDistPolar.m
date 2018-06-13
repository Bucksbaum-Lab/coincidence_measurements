function[Xp, Yp, Noutp, rs] = momentum2dDistPolar ...
            (parallel_proj, perpendicular_proj, rdata, adata, plotname, name1, name2, name3)

[Xp, Yp, Noutp, rs] = hist3polar (real(parallel_proj), real(perpendicular_proj), rdata, adata);
                
pl = pcolor(Xp, Yp, Noutp);
axis equal tight

title(plotname)
xlabel(['$$ \vec{P}(' name1 ')\parallel \left(\vec{P}(' name3 ') - \vec{P}(' name2 ')\right) $$'], 'Interpreter','latex')
ylabel(['$$ \vec{P}(' name1 ')\perp  \left(\vec{P}('    name3 ') - \vec{P}(' name2 ')\right) $$'], 'Interpreter','latex')

end