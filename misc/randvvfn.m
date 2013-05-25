function values = randvvfn(points)

path(path, './toolbox_diffc/toolbox_diffc');
path(path, './toolbox_diffc/toolbox_diffc/toolbox');
rng('default')
n = 20;
if size(points,2) == 2
    v = randn(n,n,2);
    v = perform_vf_normalization( perform_blurring(v,100) );
    [x y] = meshgrid(1:n);
    v1 = v(:,:,1);
    v2 = v(:,:,2);
    F1 = TriScatteredInterp(x(:), y(:), v1(:));
    F2 = TriScatteredInterp(x(:), y(:), v2(:));
    values = [F1(points(:,1),points(:,2)) F2(points(:,1),points(:,2))];
elseif size(points,2) == 3
    v = randn(n,n,n,3);
    v = perform_vf_normalization( perform_blurring(v,20) );
    [x y z] = meshgrid(1:n);
    v1 = v(:,:,:,1);
    v2 = v(:,:,:,2);
    v3 = v(:,:,:,3);
    F1 = TriScatteredInterp(x(:), y(:), z(:), v1(:));
    F2 = TriScatteredInterp(x(:), y(:), z(:), v2(:));
    F3 = TriScatteredInterp(x(:), y(:), z(:), v3(:));
    values = [F1(points(:,1),points(:,2),points(:,3)) F2(points(:,1),points(:,2),points(:,3)) F3(points(:,1),points(:,2),points(:,3))];
end