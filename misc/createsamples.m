samplvl = linspace(12,4,21);
x=[];
y=[];
ycurr = -10;
for i=1:numel(samplvl)
    toadd = linspace(-10,10,samplvl(i));
    y = [y repmat(ycurr,1, numel(toadd))];
    x = [x toadd];
    ycurr = ycurr+1;
end
% clearvars i samplvl toadd ycurr
points = [x' y'];
z = sin(x).*sin(y);
values = z';
