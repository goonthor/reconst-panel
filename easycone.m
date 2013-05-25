function cones = easycone(x, y, z, u, w, v, varargin)
% xmin = min(x(:));
% xmax = max(x(:));
% ymin = min(y(:));
% ymax = max(y(:));
% zmin = min(z(:));
% zmax = max(z(:));

scale = 1;
if nargin == 7
    scale = varargin{1};
end

% nconesx = 8;
% nconesy = 8;
% nconesz = 8;

% if nargin >= 9
%     nconesx = varargin{1};
%     nconesy = varargin{2};
%     nconesz = varargin{3};
% end

% if nargin == 10
%     scale = varargin{4};
% end

% xrange = linspace(xmin,xmax,nconesx);
% yrange = linspace(ymin,ymax,nconesy);
% zrange = linspace(zmin,zmax,nconesz);

% [cx cy cz] = meshgrid(xrange,yrange,zrange);

cones = coneplot(x,y,z,u,v,w,scale,'nointerp');
set(cones,'FaceColor','red','EdgeColor','none')
xlabel('X1');ylabel('X2');zlabel('X3');
view(3)
camproj perspective; 
camzoom(1.5)
camlight right; lighting phong
set(cones,'DiffuseStrength',.8) 
axis square