function bootstrpkmplot(parallel, perpendicular, nbins, title)

boot = bootstrp(1000,@(parallel,perpendicular) getcentroidbykm(parallel,perpendicular),parallel,perpendicular);

[angle1,angle1X] = hist(boot(:,4),nbins);
[angle2,angle2X] = hist(boot(:,5),nbins);
[angle3,angle3X] = hist(boot(:,6),nbins);

[R1,R1X] = hist(boot(:,1),nbins);
[R2,R2X] = hist(boot(:,2),nbins);
[R3,R3X] = hist(boot(:,3),nbins);

figure();
plot(angle1X,angle1,'-o',angle2X,angle2,'-o',angle3X,angle3,'-o')
xlabel('angle (rad)')

end
