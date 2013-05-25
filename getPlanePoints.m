function [x y z flag] = getPlanePoints(p, n, bounds, gridsize)
if n(1)==0 && n(2)==0 && n(3)==0
    disp('Error: normal cannot be all zeros.')
    flag = 0;
elseif n(1)==0 && n(2)==0
    [x y] = meshgrid(linspace(bounds(1),bounds(2),gridsize),...
                     linspace(bounds(3),bounds(4),gridsize));
    x = x(:);
    y = y(:);
    z = zeros(numel(x),1)+p(3)/n(3);
    flag = 1;
elseif n(1)==0 && n(3)==0
    [x z] = meshgrid(linspace(bounds(1),bounds(2),gridsize),...
                     linspace(bounds(5),bounds(6),gridsize));
    x = x(:);
    z = z(:);
    y = zeros(numel(x),1)+p(2)/n(2);
    flag = 2;
elseif n(2)==0 && n(3)==0
    [y z] = meshgrid(linspace(bounds(3),bounds(4),gridsize),...
                     linspace(bounds(5),bounds(6),gridsize));    
    y = y(:);
    z = z(:);
    x = zeros(numel(y),1)+p(1)/n(1);
    flag = 3;
elseif n(1)==0
    [x z] = meshgrid(linspace(bounds(1),bounds(2),gridsize),...
                     linspace(bounds(5),bounds(6),gridsize));
    x = x(:);
    z = z(:);
    y = p(2)-(n(3)/n(2))*(z-p(3));
    flag = 1;
elseif n(2)==0
    [y z] = meshgrid(linspace(bounds(3),bounds(4),gridsize),...
                     linspace(bounds(5),bounds(6),gridsize));  
    y = y(:);
    z = z(:);    
    x = p(1)-(n(3)/n(1))*(z-p(3));
    flag = 1;
elseif n(3)==0
    [x z] = meshgrid(linspace(bounds(1),bounds(2),gridsize),...
                     linspace(bounds(5),bounds(6),gridsize));
    x = x(:);
    z = z(:);    
    y = p(2)-(n(1)/n(2))*(x-p(1));
    flag = 2;
else
    [x y] = meshgrid(linspace(bounds(1),bounds(2),gridsize),...
                     linspace(bounds(3),bounds(4),gridsize));
    x = x(:);
    y = y(:);
    z = -(n(1)/n(3))*(x-p(1)) - (n(2)/n(3))*(y-p(2)) + p(3);
    flag = 1;
end