function [guesses fsize model modelfactor] = svr(points, values, newpoints, weights, modelIn, options)

pnt_dims = size(points,2);
val_dims = size(values,2);

doXVal = options.xval;
fsize = options.fsize;

if size(fsize,2) == pnt_dims
    disp('Using supplied feature sizes.')
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

[x newpnts gamma scalar] = scaledomain(points, newpoints, fsize);

if ~options.rebuildmodel
    for j=1:size(points,2)
        newpnts(:,j) = newpoints(:,j)*options.modelfactor.scalar(j);
    end
end

guesses = zeros(size(newpnts,1),val_dims);

if ~isa(modelIn, 'double')
    model = modelIn;
end

if ~isa(options.modelfactor, 'double')
    modelfactor = options.modelfactor;
end

if options.rebuildmodel
    modelfactor.scalar = scalar;
end

for i=1:val_dims
    
    % Build the model
    if options.rebuildmodel || val_dims > 1
        y = values(:, i);
        disp(['Performing SVR on dim ' num2str(i) '...']); drawnow('update')
        switch i
            case 1  
                [y modelfactor.x1] = normalize(y);
                model.x1 = svmtrain(weights, y, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1']);
            case 2
                [y modelfactor.x2] = normalize(y);
                model.x2 = svmtrain(weights, y, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1']);
            case 3
                [y modelfactor.x3] = normalize(y);
                model.x3 = svmtrain(weights, y, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1']);
        end
        
        if doXVal
            err = svmtrain(weights, y, x, ['-q -s 4 -t 2 -g ' num2str(gamma) ' -c 1 -v 5']);
            disp(['The mean squared error is ' num2str(err)])
        end
    end
    
    % Query the model
    disp('Querying model for values...'); drawnow('update')
    switch i
        case 1
            guesses(:,i) = denormalize(svmpredict(rand(size(newpnts,1),1),...
                newpnts, model.x1),modelfactor.x1);
        case 2
            guesses(:,i) = denormalize(svmpredict(rand(size(newpnts,1),1),...
                newpnts, model.x2),modelfactor.x2);
        case 3
            guesses(:,i) = denormalize(svmpredict(rand(size(newpnts,1),1),... 
                newpnts, model.x3),modelfactor.x3);
    end
end