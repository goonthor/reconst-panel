function [data factor] = normalize(d)
factor = max(abs(d(:)));
data = d./factor;