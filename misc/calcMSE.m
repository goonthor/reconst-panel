function MSE = calcMSE(x1, x2)
MSE = sum((x1-x2).^2)./numel(x1);