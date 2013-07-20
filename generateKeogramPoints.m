function newpoints = generateKeogramPoints(handles)
x1 = str2double(get(handles.editx1, 'String'));
y1 = str2double(get(handles.edity1, 'String'));
z1 = str2double(get(handles.editz1, 'String'));

x2 = str2double(get(handles.editx2, 'String'));
y2 = str2double(get(handles.edity2, 'String'));
z2 = str2double(get(handles.editz2, 'String'));

mint = min(handles.pnts(:,handles.pntdim));
maxt = max(handles.pnts(:,handles.pntdim));

N = 500; % *** ATTN: MAGIC # 500 ***

t = linspace(mint,maxt,N);
t = repmat(t,N,1);
t = t(:);

if handles.pntdim == 3
    [x y] = getNPointsFromSegment2D(x1,y1,x2,y2,N);
    x = repmat(x,N,1);
    y = repmat(y,N,1);
    newpoints = [x y t];
else
    [x y z] = getNPointsFromSegment3D(x1,y1,z1,x2,y2,z2,N);
    x = repmat(x,N,1);
    y = repmat(y,N,1);
    z = repmat(z,N,1);
    newpoints = [x y z t];
end

function [x y] = getNPointsFromSegment2D(x1,y1,x2,y2,N)
t = linspace(0,1,N)';
x = x1+(x2-x1)*t;
y = y1+(y2-y1)*t;

function [x y z] = getNPointsFromSegment3D(x1,y1,z1,x2,y2,z2,N)
t = linspace(0,1,N)';
x = x1+(x2-x1)*t;
y = y1+(y2-y1)*t;
z = z1+(z2-z1)*t;