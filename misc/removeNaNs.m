function [points_p values_p] = removeNaNs(points, values)
[rp ~] = find(isnan(points));
[rv ~] = find(isnan(values));
r = [rp; rv];
r = unique(r);
points_p = points;
values_p = values;
points_p(r,:) = [];
values_p(r,:) = [];
disp(['Removed ' num2str(numel(r)) ' bad points.'])