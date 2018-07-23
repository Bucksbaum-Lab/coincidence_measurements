function RanglePair = getcentroidbykm(parallel, perpendicular) 

[~, centroids] = kmeans([parallel,perpendicular],3);

angle = atan2(centroids(:,2),centroids(:,1));
[angle,I] = sort(angle);
R = sqrt(centroids(:,1).^2+centroids(:,2).^2);
R = R(I,:);

RanglePair = [R,angle];

end