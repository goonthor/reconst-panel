function [guesses fsize] = svr_window(points, values, newpoints, weights, options)

pnt_dims = size(points,2);
val_dims = size(values,2);

doXVal = options.xval;
fsize = options.fsize;

if size(fsize,2) == pnt_dims
    disp('Using supplied feature sizes: ')
else
    disp('Using default feature size.')
    fsize = zeros(1, pnt_dims);
    fsize = fsize + (max(points(:))-min(points(:)))/10;
end

if numel(find(fsize==0)) > 0
    disp('Feature size cannot be zero. Using default feature size.')
    fsize = zeros(1, pnt_dims);
    fsize = fsize + (max(points(:))-min(points(:)))/10;
end

if numel(weights) ~= numel(points(:,1))
    disp('Using equal weights- no valid weights file specified.')
    weights = ones(numel(points(:,1)),1);
end

window = fsize(pnt_dims);

guesses = zeros(size(newpoints,1),val_dims);
for i=1:val_dims
    y = values(:, i);
    fprintf(['Performing windowed SVR on dim ' num2str(i) '.'])
    guesses(:,i) = windowhelper(points, y, newpoints, weights, fsize, doXVal, window);
end

function dim_guesses = windowhelper(points, y, newpoints, weights, fsize, doXVal, window)
dim_guesses = zeros(size(newpoints,1),1);
pntdim = size(points,2);
acc=0;
curr = min(newpoints(:,pntdim));
times = 0;
badtimes = [];
while numel(curr) ~= 0
    currind = newpoints(:,pntdim) == curr;
    pnt_ind = find(points(:,pntdim) >= curr-window & points(:,pntdim) <= curr+window);
    if numel(pnt_ind)>0
        [y_norm factor] = normalize(y(pnt_ind));
        
        [x newpnts gamma] = scaledomain(points(pnt_ind,:), newpoints(currind,:), fsize);
        
%         disp(['  Building model for time ' num2str(curr) ' from ' num2str(points(min(pnt_ind),pntdim)) ' to ' num2str(points(max(pnt_ind),pntdim)) ' (' num2str(numel(pnt_ind)) ').']); drawnow('update')
        model = svmtrain(weights(pnt_ind),y_norm, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1']);
        if doXVal
            acc = acc + svmtrain(weights(pnt_ind),y_norm, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1 -v 5']);
        end
%         disp(['  Predicting ' num2str(sum(currind)) ' values for model.']); drawnow('update')
        [temp] = svmpredict(rand(size(newpnts,1),1), newpnts, model);
        dim_guesses(currind) = denormalize(temp,factor);
        times = times + 1;
    else
        badtimes = [badtimes curr];
    end
    if mod(times,20)==0
        disp('.')
    else
        fprintf('.')
    end
    ind = newpoints(:,pntdim) > curr;
    curr = min(newpoints(ind,pntdim));
    
end
disp('.')
disp(['Number of analyses done: ' num2str(times)])

if numel(badtimes) > 0
    disp('Models could not be produced for the following times:')
    disp(badtimes)
    disp('To rectify this, increase the time feature size.')
end

if doXVal
    disp(['Average cross-validation error: ' num2str(acc/times)])
end
