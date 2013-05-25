function [guesses fsize] = svr(points, values, new_points, varargin)

pnt_dims = size(points,2);
val_dims = size(values,2);

fsize = zeros(1, pnt_dims);
fsize = fsize + (max(points(:))-min(points(:)))/10;

split=1;

i=1;
while i<=nargin-3
    switch varargin{i}
        case 'fsize'
            if size(varargin{i+1},2) == pnt_dims
                fsize = varargin{i+1};
                disp('Using supplied feature sizes: ')
                disp(fsize);
            else
                disp('Using default feature size.')
            end
            i=i+1;
        case 'split'
            split = varargin{i+1};
    end
    i=i+1;
end

if numel(find(fsize==0)) > 0
    disp('Feature size cannot be zero. Using default feature size.')
    fsize = zeros(1, pnt_dims);
    fsize = fsize + (max(points(:))-min(points(:)))/10;
end

% total_domain = scaleinput([points; new_points], fsize);
% x = total_domain(1:num_samples, :);
% new_points = total_domain(num_samples+1:end, :);

guesses = zeros(size(new_points,1),val_dims);
disp(['Splitting the samples into ' num2str(split) ' models.']); drawnow('update')
for i=1:val_dims
    y = values(:, i);
    % [y factor] = normalize(y);
    disp(['Training SVM on dim ' num2str(i) '...']); drawnow('update')
    % model = svmtrain(y, x, '-q -s 4 -t 2 -g 1 -c 10');
    % disp('Predicting values...')
    % drawnow('update')
    guesses(:,i) = splithelper(points, y, new_points, split, fsize);
end

% assume: sample points sorted last dimension (the one to split)
function guesses = splithelper(points, y, new_points, split, fsize)
numsamples = numel(y);
numpermodel = floor(numsamples/split);
guesses = zeros(size(new_points,1),1);
pntdim = size(points,2);
disp(['Using ' num2str(numpermodel) ' samples per model.']); drawnow('update')
for j=1:split
    if j ~= split
        ind=(j-1)*numpermodel+1:j*numpermodel;
    else
        ind=(j-1)*numpermodel+1:numsamples;
    end
    pmin = min(points(ind,pntdim));
    pmax = max(points(ind,pntdim));
%     disp(num2str(pmin))
%     disp(num2str(pmax))
    ind2 = find(new_points(:,pntdim) >= pmin & new_points(:,pntdim) <= pmax);
%     disp(num2str(numel(new_points(:,pntdim))))
    [y_norm factor] = normalize(y(ind));

    model_domain = scaleinput([points(ind,:); new_points(ind2,:)], fsize);
    x = model_domain(1:numel(ind), :);
    newpnts = model_domain(numel(ind)+1:end, :);

    disp(['  Building model ' num2str(j) ' from samples ' num2str(min(ind)) ' to ' num2str(max(ind)) '.']); drawnow('update')
    model = svmtrain(y_norm, x, '-q -s 4 -t 2 -g 1 -c 10');
    disp(['  Predicting ' num2str(numel(ind2)) ' values for model ' num2str(j) '.']); drawnow('update')
    guesses(ind2) = denormalize(svmpredict(rand(size(newpnts,1),1), newpnts, model),factor);
end


% function guesses = splithelper(x, y, new_points, split)
% numsamples = numel(y);
% numpermodel = floor(numsamples/split);
% guesses = zeros(size(new_points,1),1);
% pntdim = size(x,2);
% disp(['  Using ' num2str(numpermodel) ' samples per model.']); drawnow('update')
% for j=1:split
%     if j ~= split
%         ind=(j-1)*numpermodel+1:j*numpermodel;
%     else
%         ind=(j-1)*numpermodel+1:numsamples;
%     end
%     disp(['  Building model ' num2str(j) '.']); drawnow('update')
%     model = svmtrain(y(ind), x(ind,:), '-q -s 4 -t 2 -g 1 -c 10');
%     disp(['  Predicting values for model ' num2str(j) '.']); drawnow('update')
%     xmin = min(x(ind,pntdim));
%     xmax = max(x(ind,pntdim));
%     ind2 = find(new_points(:,pntdim) <= xmin & new_points(:,pntdim) >= xmax);
%     guesses(ind2) = svmpredict(rand(size(new_points(ind2),1),1), new_points(ind2), model);
% end

