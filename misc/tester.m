path(path, './toolbox_diffc/toolbox_diffc');
path(path, './toolbox_diffc/toolbox_diffc/toolbox');
rng('default')
n = 20;
V = randn(n,n,n,3);
V = perform_vf_normalization( perform_blurring(V,20) );
u = V(:,:,:,1);
v = V(:,:,:,2);
w = V(:,:,:,3);
[x y z] = meshgrid(1:n);
xmin = min(x(:));
xmax = max(x(:));
ymin = min(y(:));
ymax = max(y(:));
zmin = min(z(:));
zmax = max(z(:));
xrange = linspace(xmin,xmax,4);
yrange = linspace(ymin,ymax,4);
zrange = linspace(zmin,zmax,4);
[cx cy cz] = meshgrid(xrange,yrange,zrange);

h = streamline(x,y,z,u,v,w,cx,cy,cz);
set(h,'Color','red')
view(3)
% cones = coneplot(x,y,z,u,v,w,cx,cy,cz);
% set(cones,'FaceColor','red','EdgeColor','none')
% xlabel('X1');ylabel('X2');zlabel('X3');
% view(3)
% camproj perspective;
% camzoom(1.5)
% camlight right; lighting phong
% set(cones,'DiffuseStrength',.8)
% axis square