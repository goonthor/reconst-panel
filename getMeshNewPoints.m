function pntmatrix = getMeshNewPoints(pntdim, newpnts)                    %
pntmatrix = newpnts.x1;                                                   %
if pntdim == 2                                                            %
    [X1 X2] = meshgrid(newpnts.x1, newpnts.x2);
    pntmatrix = [X1(:) X2(:)];
elseif pntdim == 3
    [X1 X2 X3] = meshgrid(newpnts.x1, newpnts.x2, newpnts.x3);
    pntmatrix = [X1(:) X2(:) X3(:)];
elseif pntdim == 4
    [X1 X2 X3 X4] = ndgrid(newpnts.x1,newpnts.x2,newpnts.x3,newpnts.x4);
    X1 = permute(X1,[2 1 3 4]);
    X2 = permute(X2,[2 1 3 4]);
    X3 = permute(X3,[2 1 3 4]);
    X4 = permute(X4,[2 1 3 4]);
    pntmatrix = [X1(:) X2(:) X3(:) X4(:)];
end