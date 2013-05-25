function h = viewPlane(x, y, z, c, flag)
switch flag
    case 1
        gridDelaunay = delaunay(x,y);
    case 2
        gridDelaunay = delaunay(x,z);
    case 3
        gridDelaunay = delaunay(y,z);
    otherwise
        disp('Error is viewing plane.')
end
h = trisurf(gridDelaunay,x,y,z,c); 
set(h,'EdgeColor','none')