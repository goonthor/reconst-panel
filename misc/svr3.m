function [guesses fsize] = svr3(points, values, newpoints, options)

pnt_dims = size(points,2);
val_dims = size(values,2);

doXVal = options.xval;
fsize = options.fsize;
split=options.split;

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

guesses = zeros(size(newpoints,1),val_dims);
for i=1:val_dims
    y = values(:, i);
    disp(['Performing split-SVR on dim ' num2str(i) '...']); drawnow('update')
    guesses(:,i) = splithelper(points, y, newpoints, split, fsize, doXVal);
end


% assume: sample points sorted last dimension (the one to split)
function dim_guesses = splithelper(points, y, new_points, split, fsize, doXVal)
numsamples = numel(y);
numpermodel = floor(numsamples/split);
dim_guesses = zeros(size(new_points,1),1);
pntdim = size(points,2);
disp(['Using ' num2str(numpermodel) ' samples per model.']); drawnow('update')
acc=0;
for j=1:split
    if j ~= split
        samp_ind=(j-1)*numpermodel+1:j*numpermodel;
    else
        samp_ind=(j-1)*numpermodel+1:numsamples;
    end
    pmin = min(points(samp_ind,pntdim));
    pmax = max(points(samp_ind,pntdim));
    new_ind = find(single(new_points(:,pntdim)) >= single(pmin) &...
        single(new_points(:,pntdim)) <= single(pmax));
    [y_norm factor] = normalize(y(samp_ind));
    
    [x newpnts gamma] = scaledomain(points(samp_ind,:), new_points(new_ind,:), fsize);

    disp(['  Building model ' num2str(j) ' from samples ' num2str(min(samp_ind)) ' to ' num2str(max(samp_ind)) '.']); drawnow('update')
    
    model = svmtrain(y_norm, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1']);
    if doXVal
        acc = acc + svmtrain(y_norm, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1 -v 5']);
    end
    disp(['  Predicting ' num2str(numel(new_ind)) ' values for model ' num2str(j) '.']); drawnow('update')
    dim_guesses(new_ind) = denormalize(svmpredict(rand(size(newpnts,1),1), newpnts, model),factor);
end
disp(['Average cross-validation accuracy is ' num2str(acc/split) '.'])
