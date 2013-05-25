
[X1 X2 X3 X4] = ndgrid(1:10,1:10,1:10,1:10);
X1 = permute(X1,[2 1 3 4]);
X2 = permute(X2,[2 1 3 4]);
X3 = permute(X3,[2 1 3 4]);
X4 = permute(X4,[2 1 3 4]);
points = [X1(:) X2(:) X3(:) X4(:)];
values = zeros(10000,1);
ind = find((points(:,1)==5 | points(:,1)==6) &...
    (points(:,2)==5 | points(:,2)==6) &...
    (points(:,3)==5 | points(:,3)==6) &...
    (points(:,4)==5 | points(:,4)==6));
values(ind) = 1;