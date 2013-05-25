function [u v w f] = vv4dfn(x,y,z,t)
f = x.*sin(y.*z+t/10);
u = sin(y.*z+t/10);
v = x.*z.*cos(y.*z+t/10);
w = x.*y.*cos(y.*z+t/10);