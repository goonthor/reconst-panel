function [newpoints flag] = generateSlicePoints(handles)
x1 = str2double(get(handles.editx1, 'String'));
y1 = str2double(get(handles.edity1, 'String'));
z1 = str2double(get(handles.editz1, 'String'));

x2 = str2double(get(handles.editx2, 'String'));
y2 = str2double(get(handles.edity2, 'String'));
z2 = str2double(get(handles.editz2, 'String'));

mint = min(handles.pnts(:,handles.pntdim));
maxt = max(handles.pnts(:,handles.pntdim));

Nt = 30;
Ns = 50;
t = linspace(mint,maxt,Nt);
% need to increase Nt?

t = repmat(t,Ns^2,1);
t = t(:);
bounds = [min(handles.pnts(:,1)) max(handles.pnts(:,1)) min(handles.pnts(:,2)) max(handles.pnts(:,2)) min(handles.pnts(:,3)) max(handles.pnts(:,3))];
[x y z flag] = getPlanePoints([x1 y1 z1], [x2 y2 z2], bounds, Ns);
x = repmat(x,Nt,1);
y = repmat(y,Nt,1);
z = repmat(z,Nt,1);

newpoints = [x y z t];

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