function [p_sc newp_sc gamma scalar] = scaledomain(points, newpoints, fsize)
% 1. obtain factor to scale each dimension to [-1, 1]
% disp('factor=')
factor = max(max(abs(points)),max(abs(newpoints)));
% factor = [6.282722513089006 72 .4039];

% 2. scale the feature sizes as well
fsize = fsize./factor;

% 3. the fwhm used to calculate gamma is the max of these because features
% sizes are proxies for FWHM's.  The largest FWHM relative to the range of
% the domain calculates gamma. The other dimensions are scaled to greater
% lengths.
fwhm = max(fsize);

% 4. The scale lengths range from 1 to inf where 1 is the dimension with
% the FWHM. inf signifies very small fsize relative to domain length so
% this dimension must be stretched more.
scalar = fwhm./fsize;

% ** The initial scaling factor (to [-1, 1]) is built in here
scalar = scalar./factor;

% Gamma is parameter of gaussian = e^(-gamma*|u-v|^2)
gamma = 2*log(2)/(fwhm^2);

p_sc = zeros(size(points));
newp_sc = zeros(size(newpoints));

% Scale here.
for j=1:size(points,2)
    p_sc(:,j) = points(:,j)*scalar(j);
    newp_sc(:,j) = newpoints(:,j)*scalar(j);
end


% disp('Scalar vector is: ')
% scalar
% disp(['FWHM is ' num2str(fwhm)])
% disp(['Gamma is ' num2str(gamma)])
