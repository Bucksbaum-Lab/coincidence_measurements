function RanglePair = getcentroidbyhir(parallel, perpendicular)

T = clusterdata([parallel,perpendicular],'linkage','average','Maxclust',3);

centroid1 = [mean(parallel(T==1)),mean(perpendicular(T==1))];
centroid2 = [mean(parallel(T==2)),mean(perpendicular(T==2))];
centroid3 = [mean(parallel(T==3)),mean(perpendicular(T==3))];

centroids = [centroid1;centroid2;centroid3];

angle = atan2(centroids(:,2),centroids(:,1));
[angle,I] = sort(angle);
R = sqrt(centroids(:,1).^2+centroids(:,2).^2);
R = R(I,:);

RanglePair = [R,angle];

end