function easyslice(x, y, z, v)
rx = max(x(:))-min(x(:));
ry = max(y(:))-min(y(:));
rz = max(z(:))-min(z(:));
xslice = [min(x(:))+rx/20, min(x(:))+rx/3, min(x(:))+2*rx/3, max(x(:))-rx/20];
yslice = max(y(:))-ry/20;
zslice = [min(z(:))+rz/20, min(z(:))+rz/3, min(z(:))+2*rz/3,];
slice(x,y,z,v, xslice, yslice, zslice);
axis square
grid off
set(findobj(gca,'Type','Surface'),'EdgeColor','none')
xlabel('X1');ylabel('X2');zlabel('X3');
alpha(.75)