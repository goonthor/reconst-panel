function visualize(points, values, newpoints, guesses, varargin)

% Set default values for these variables.
isAnimated = 0;
showSamples = 0;
showMag = 0;
arrowsize = 1;
isVector = size(values,2)-1;
currp = 1;
fps = 2;
sbs = 0;
sampID = 1;
dispID = 1;
titlep1 = 'No Title';
cbounds = [-1 1];
isKeo = 0;
isSlice = 0;

% Parse the input
for i=1:nargin-4
    switch varargin{i}
        case 'animate'
            isAnimated = varargin{i+1};
        case 'showSamples'
            showSamples = varargin{i+1};
        case 'showMag'
            showMag = varargin{i+1};
        case 'arrowsize'
            arrowsize = varargin{i+1};
        case 'panel'
            currp = varargin{i+1};
        case 'fps'
            fps = varargin{i+1};
        case 'sbs'
            sbs = varargin{i+1};
        case 'dispID'
            dispID = varargin{i+1};
        case 'sampID'
            sampID = varargin{i+1};
        case 'title'
            titlep1 = varargin{i+1};
        case 'keo'
            isKeo = varargin{i+1};
        case 'slice'
            isSlice = varargin{i+1};
        case 'cbounds'
            cbounds(1) = varargin{i+1};
            cbounds(2) = varargin{i+2};
    end
end

% Make sure arrow size is valid.
if isnan(arrowsize)
    arrowsize=1;
end

% Bring up the window and get the current view orientation.
figure(1)
[az el] = view;

pntdim = size(points, 2);
numpnts = size(points, 1);

% This function is used to switch the coordinate system of the points and
% newpoints.
switcherfn = @(X1,X2,X3) switchcoordsys(dispID,sampID,X1,X2,X3,pntdim);


% This is the beginning of the main block of code. There are 4 possible
% branches of this "if" statement: keogram, slice plane, mesh data, and
% scattered data.
% *** This is the Keogram branch.
if isKeo
    figure(2)
    clf(2)
    
    % Switch coordinate systems if necessary.
    [n1 n2 n3]=switcherfn(newpoints(:,1),newpoints(:,2),newpoints(:,3));
    newpoints(:,1) = n1;
    newpoints(:,2) = n2;
    newpoints(:,3) = n3;
    
    % ATTN: MAGIC # 500 below
    % SEE generateKeogramPoints
    [X1 X2] = meshgrid(unique(newpoints(:,pntdim)), 1:500); % ***
    g = reshape(guesses,size(X1));
    h = surf(X1,X2,g);
    set(h,'EdgeColor','none')
    xlabel('Time')
    ylabel('Keogram Line');
    view([az el])
    colorbar;
% *** This is the slice plane branch.
elseif isSlice
    figure(2)
    clf(2)
    
    % TODO: add 3D option here.
    % if pntdim == 4
    % Switch the coordinate systems if necessary.
    [n1 n2 n3]=switcherfn(newpoints(:,1),newpoints(:,2),newpoints(:,3));
    newpoints(:,1) = n1;
    newpoints(:,2) = n2;
    newpoints(:,3) = n3;
    [p1 p2 p3]=switcherfn(points(:,1),points(:,2),points(:,3));
    points(:,1)=p1;
    points(:,2)=p2;
    points(:,3)=p3;
    
    % Sorted list of all times.
    times = sort(unique(newpoints(:,4)));
    
    % Frames in the animation.
    numframes = numel(times);
    
    % This varible determines which samples are shown in the current frame.
    % Assuming the spacing between times is constant (which we know to be
    % true since we generate them in generateSlicePoints.m) we can take the
    % difference between the first two times and divide it by two.
    slack = (times(2)-times(1))/2;
    
    % Number of new points in a frame.
    pntsframe = numel(newpoints(:,1))/numframes;
    
    % Find the bounds so the axes and color scale does not change at each
    % frame of the animation.
    xmin = min(points(:,1));
    xmax = max(points(:,1));
    ymin = min(points(:,2));
    ymax = max(points(:,2));
    zmin = min(points(:,3));
    zmax = max(points(:,3));
    bounds = [xmin xmax ymin ymax zmin zmax cbounds(1) cbounds(2)];
    
    set(gca,'NextPlot','replaceChildren');
    
    for j=1:numframes
        title({titlep1; ['Time: ' num2str(times(j))]})
        
        % Current index in the list of new points.
        curr = (j-1)*pntsframe+1;
        
        viewPlane(newpoints(curr:curr+pntsframe-1,1),...
            newpoints(curr:curr+pntsframe-1,2),...
            newpoints(curr:curr+pntsframe-1,3),...
            guesses(curr:curr+pntsframe-1), isSlice);
        
        % Find any samples which are +/- slack distance from the current
        % time in the animation and show them.
        if showSamples
            ind = find(points(:,4) > times(j)-slack &... 
                points(:,4) < times(j)+slack);
            hold on
            scatter3(points(ind,1),points(ind,2), points(ind,3),[],...
                values(ind),'filled');
            hold off
        end
        xlabel('X1')
        ylabel('X2')
        zlabel('X3')
        axis(bounds);
        view([az el])
        pause(1/fps);
    end
    % else
        % Write the code for 3D slice plane here.
    % end
% This is the mesh branch. In this branch and the scattered branch, we
% switch on the dimensionality of the data.  These branches could probably
% be refactored to become a little cleaner, but there are enough minor
% differences between them that refactoring is difficult.
% Within each dimension we have to test if the data is vector or scalar and
% if it needs to be animated or not.  Vector animations have not been
% implemented.
elseif isstruct(newpoints)
    switch pntdim
        % 1-D data.
        case 1
            figure(1)
            clf(1)
            hold on
            plot(newpoints.x1, guesses)
            if showSamples
                plot(points, values, 'r*')
            end
            title({titlep1; ['Reconstructed Signal with ',...
                num2str(numpnts), ' Samples']});
            hold off
        % 1-D data.
        case 2
            % Create the mesh.
            [X1 X2] = meshgrid(newpoints.x1, newpoints.x2);
            
            figure(1)
            clf(1)
            
            % Switch coordinate systems if necessary.
            [X1 X2] = switcherfn(X1, X2, 0);
            [p1 p2] = switcherfn(points(:,1),points(:,2), 0);
            points(:,1)=p1;
            points(:,2)=p2;
            
            % Vector data.
            if isVector
                hold on
                quiver(X1(:), X2(:), guesses(:,1), guesses(:,2), arrowsize);
                if showMag
                    mag = sqrt(guesses(:,1).^2 + guesses(:,2).^2);
                    h=surf(X1, X2, reshape(mag,size(X1)));
                    set(h,'EdgeColor','none')
                    alpha(.5)
                    % colorbar;
                end
                if showSamples
                    quiver(points(:,1),points(:,2), values(:,1), values(:,2), .5, 'color', [1 0 0]);
                end
                title({titlep1; ['Reconstructed Signal With ', num2str(numpnts), ' Samples']});
                xlabel('X1')
                ylabel('X2')
                hold off
            else
                % Reshape the guesses vector so meesh visualizations can be
                % used, i.e. surf.
                g = reshape(guesses,size(sX1));
                
                % Animation
                if isAnimated
                    
                    % Find the bounds so the axes and color scale does not change at each
                    % frame of the animation.
                    xmin = min(X1(:));
                    xmax = max(X1(:));
                    ymin = min(guesses);
                    ymax = max(guesses);
                    bounds = [xmin xmax ymin ymax];

                    set(gca,'NextPlot','replaceChildren');
                    numframes = numel(newpoints.x2);
                    for j=1:numframes
                        title({titlep1; ['Time: ' num2str(newpoints.x2(j))]})
                        plot(X1(j,:),g(j,:));
                        axis(bounds);
                        pause(1/fps);
                    end
                % Frame-by-frame
                elseif sbs
                    plottingfn = @(ind) plothelper(X1,g,ind);
                    panelplotter(numel(newpoints.x2),currp,plottingfn,newpoints.x2,az,el);
                % 2D surface.
                else
                    hold on
                    h=surf(X1,X2,g);
                    set(h,'EdgeColor','none')
                    if showSamples
                        plot3(points(:,1),points(:,2), values, 'r*');
                    end
                    title({titlep1; ['Reconstructed Signal With ', num2str(numpnts), ' Samples']});
                    colorbar;
                    xlabel('X1')
                    ylabel('X2')
                    hold off
                end
            end
        % 3-D data.
        case 3
            % Create the mesh.
            [X1 X2 X3] = meshgrid(newpoints.x1, newpoints.x2, newpoints.x3);
            
            figure(1)
            clf(1)
            
            % Switch coordinate systems if necessary.
            [X1 X2 X3]=switcherfn(X1,X2,X3);
            [p1 p2 p3]=switcherfn(points(:,1),points(:,2),points(:,3));
            points(:,1)=p1;
            points(:,2)=p2;
            points(:,3)=p3;
            
            % Vector data
            if isVector
                hold on
                easycone(X1(:), X2(:), X3(:),...
                    guesses(:,1), guesses(:,2), guesses(:,3),...
                    arrowsize);
                if showSamples
                    quiver3(points(:,1),points(:,2), points(:,3), values(:,1), values(:,2), values(:,3));
                end
                if showMag
                    mag = sqrt(guesses(:,1).^2 + guesses(:,2).^2 + guesses(:,3).^2);
                    easyslice(X1, X2, X3, reshape(mag,size(X1)));
                end
                title({titlep1; ['Reconstructed Signal With ', num2str(numpnts), ' Samples']});
                xlabel('X1')
                ylabel('X2')
                zlabel('X3')
                hold off
            else
                % Reshape the guesses vector so meesh visualizations can be
                % used, i.e. surf.
                g = reshape(guesses,size(X1));
                
                % Animation
                if isAnimated
                    % Find the bounds so the axes and color scale does not change at each
                    % frame of the animation.
                    xmin = min(X1(:));
                    xmax = max(X1(:));
                    ymin = min(X2(:));
                    ymax = max(X2(:));
                    zmin = min(guesses);
                    zmax = max(guesses);
                    bounds = [xmin xmax ymin ymax zmin zmax cbounds(1) cbounds(2)];

                    set(gca,'NextPlot','replaceChildren');
                    numframes = numel(newpoints.x3);
                    
                    % This varible determines which samples are shown in 
                    % the current frame. Assuming the spacing between times 
                    % is constant (which we know to be true since we 
                    % generate them in generateSlicePoints.m) we can take 
                    % the difference between the first two times and divide 
                    % it by two.
                    slack = (newpoints.x3(2)-newpoints.x3(1))/2;
                    for j=1:numframes
                        title({titlep1; ['Time: ' num2str(newpoints.x3(j))]})
                        h = surf(X1(:,:,j),X2(:,:,j),g(:,:,j));
                        set(h,'EdgeColor','none')
                        if showSamples
                            ind = find(points(:,3) > newpoints.x3(j)-slack & points(:,3) < newpoints.x3(j)+slack);
                            hold on
                            % Colored points.
                            % scatter3(points(ind,1),points(ind,2), values(ind), 15,values(ind), 'filled');
                            plot3(points(ind,1),points(ind,2), values(ind), 'r*');
                            hold off
                        end
                        view([az el])
                        axis(bounds);
                        pause(1/fps);
                    end
                elseif sbs
                    if showSamples
                        disp('Due to computational constraints, samples cannot be shown.')
                    end
                    plottingfn = @(ind) surfhelper(X1,X2,g,cbounds,ind);
                    panelplotter(numel(newpoints.x3),currp,plottingfn,newpoints.x3,az,el);
                else
                    hold on
                    if dispID==1 && sampID==1
                        easyslice(X1, X2, X3, g);
                    else
                        scatter3(X1(:), X2(:), X3(:), 4, guesses,'filled');
                    end
                    if showSamples
                        scatter3(points(:,1),points(:,2), points(:,3), 15,values, 'filled');
                    end
                    title({titlep1; ['Reconstructed Signal With ',num2str(numpnts),' Samples']});
                    colorbar();
                    xlabel('X1')
                    ylabel('X2')
                    zlabel('X3')
                    hold off
                end
            end
        case 4
            figure(1)
            clf(1)
            % meshgrid only works up to 3D, so use ndgrid instead.
            [X1 X2 X3 ~] = ndgrid(newpoints.x1, newpoints.x2, newpoints.x3, newpoints.x4);
            
            % The formatting with ndgrid is slightly different, so we have
            % to permute the matrices.
            X1 = permute(X1,[2 1 3 4]);
            X2 = permute(X2,[2 1 3 4]);
            X3 = permute(X3,[2 1 3 4]);
            
            % Switch coordinate systems if necessary.
            [X1 X2 X3]=switcherfn(X1,X2,X3);
            if showSamples
                disp('Due to computational constraints, samples cannot be shown.')
            end
            if isVector
                if isAnimated
                    set(gca,'NextPlot','replaceChildren');
                    numframes = numel(newpoints.x4);
                    for j=1:numframes
                        title({titlep1; ['Time: ' num2str(newpoints.x4(j))]})
                        easycone(X1(:), X2(:), X3(:),...
                            guesses(:,1), guesses(:,2), guesses(:,3),...
                            arrowsize);
                        pause(1/fps);
                    end
                else
                    g1 = reshape(guesses(:,1),size(X1));
                    g2 = reshape(guesses(:,2),size(X1));
                    g3 = reshape(guesses(:,3),size(X1));
                    plottingfn = @(ind) easyconehelper(X1,X2,X3,g1,g2,g3,arrowsize,ind);
                    panelplotter(numel(newpoints.x4),currp,plottingfn,newpoints.x4,az,el);
                end
            else
                guesses_reshp = reshape(guesses,size(X1));
                if isAnimated
                    % Find the bounds so the axes and color scale does not change at each
                    % frame of the animation.
                    xmin = min(X1(:));
                    xmax = max(X1(:));
                    ymin = min(X2(:));
                    ymax = max(X2(:));
                    zmin = min(X3(:));
                    zmax = max(X3(:));
                    bounds = [xmin xmax ymin ymax zmin zmax cbounds(1) cbounds(2)];
                    
                    set(gca,'NextPlot','replaceChildren');
                    numframes = numel(newpoints.x4);
                    for j=1:numframes
                        easyslice(X1(:,:,:,j), X2(:,:,:,j), X3(:,:,:,j), guesses_reshp(:,:,:,j));
                        axis(bounds)
                        title({titlep1; ['Time: ' num2str(newpoints.x4(j))]})
                        view([az el])
                        pause(1/fps)
                    end
                else
                    plottingfn = @(ind) easyslicehelper(X1,X2,X3,guesses_reshp,ind);
                    panelplotter(numel(newpoints.x4),currp,plottingfn,newpoints.x4,az,el);
                end
            end
            
        otherwise
            disp('Dimension is too great to visualize');
    end
% This is the scattered branch. In this branch and the mesh branch, we
% switch on the dimensionality of the data.  These branches could probably
% be refactored to become a little cleaner, but there are enough minor
% differences between them that refactoring is difficult.
% Within each dimension we have to test if the data is vector or scalar and
% if it needs to be animated or not.  Vector animations have not been
% implemented.
else
    switch pntdim
        case 1
            figure(1)
            clf(1)
            hold on
            scatter(newpoints, guesses)
            if showSamples
                plot(points, values, 'r*')
            end
            title({titlep1; ['Reconstructed Signal with ', num2str(numpnts), ' Samples']});
            hold off
        % Need to implement animation for this case still.
        case 2
            figure(1)
            clf(1)
            % Switch the coordinate systems if necessary.
            [n1 n2]=switcherfn(newpoints(:,1), newpoints(:,2),0);
            newpoints(:,1) = n1;
            newpoints(:,2) = n2;
            [p1 p2]=switcherfn(points(:,1),points(:,2),0);
            points(:,1)=p1;
            points(:,2)=p2;
            
            if isVector
                hold on
                quiver(newpoints(:,1), newpoints(:,2), guesses(:,1), guesses(:,2), arrowsize);
                if showMag
                    mag = sqrt(guesses(:,1).^2 + guesses(:,2).^2);
                    gridDelaunay = delaunay(newpoints(:,1), newpoints(:,2));
                    h = trisurf(gridDelaunay, newpoints(:,1), newpoints(:,2), mag);
                    set(h,'EdgeColor','none')
                end
                if showSamples
                    quiver(points(:,1),points(:,2), values(:,1), values(:,2), 0, 'color', [1 0 0]);
                end
                title({titlep1; ['Reconstructed Signal With ', num2str(numpnts), ' Samples']});
                hold off
            else
                hold on
                gridDelaunay = delaunay(newpoints(:,1), newpoints(:,2));
                h = trisurf(gridDelaunay, newpoints(:,1), newpoints(:,2), guesses);
                set(h,'EdgeColor','none')
                if showSamples
                    plot3(points(:,1),points(:,2), values, 'r*');
                end
                title({titlep1; ['Reconstructed Signal With ', num2str(numpnts), ' Samples']});
                colorbar();
                hold off
            end
        case 3
            figure(1)
            clf(1)
            % Switch the coordinate systems if necessary.
            [n1 n2 n3]=switcherfn(newpoints(:,1),newpoints(:,2),newpoints(:,3));
            newpoints(:,1) = n1;
            newpoints(:,2) = n2;
            newpoints(:,3) = n3;
            [p1 p2 p3]=switcherfn(points(:,1),points(:,2),points(:,3));
            points(:,1)=p1;
            points(:,2)=p2;
            points(:,3)=p3;
            
            if isVector
                hold on
                easycone(newpoints(:,1), newpoints(:,2), newpoints(:,3),...
                    guesses(:,1), guesses(:,2), guesses(:,3),...
                    arrowsize);
                if showSamples
                    quiver3(points(:,1),points(:,2), points(:,3), values(:,1), values(:,2), values(:,3));
                end
                if showMag
                    mag = sqrt(guesses(:,1).^2 + guesses(:,2).^2 + guesses(:,3).^2);
                    scatter3(newpoints(:,1), newpoints(:,2), newpoints(:,3), 15, mag, '*');
                end
                view([az el])
                title({titlep1; ['Reconstructed Signal With ', num2str(numpnts), ' Samples']});
                xlabel('X1')
                ylabel('X2')
                zlabel('X3')
                hold off
            else
                if isAnimated
                    % Find the bounds so the axes and color scale does not change at each
                    % frame of the animation.
                    xmin = min(newpoints(:,1));
                    xmax = max(newpoints(:,1));
                    ymin = min(newpoints(:,2));
                    ymax = max(newpoints(:,2));
                    zmin = min(guesses);
                    zmax = max(guesses);
                    bounds = [xmin xmax ymin ymax zmin zmax cbounds(1) cbounds(2)];
                    
                    set(gca,'NextPlot','replaceChildren');
                    
                    % Sorted list of all times.
                    times = sort(unique(newpoints(:,3)));
                    numframes = numel(times);
                    for j=1:numframes
                        % Find all the new points in the given time.
                        ind = find(newpoints(:,3)==times(j));
                        title({titlep1; ['Time: ' num2str(times(j))]})
                        gridDelaunay = delaunay(double(newpoints(ind,1)), double(newpoints(ind,2)));
                        h = trisurf(gridDelaunay, newpoints(ind,1), newpoints(ind,2), guesses(ind));
                        set(h,'EdgeColor','none')
                        % scatter(newpoints(ind,1), newpoints(ind,2), 4, guesses(ind),'filled');
                        view([az el])
                        axis(bounds);
                        pause(1/fps);
                    end
                elseif sbs
                    times = sort(unique(newpoints(:,3)));
                    plottingfn = @(tind) trisurfhelper(newpoints,guesses,times,tind);
                    panelplotter(numel(times),currp,plottingfn,times,az,el);
                else
                    hold on
                    scatter3(newpoints(:,1), newpoints(:,2), newpoints(:,3), 4, guesses,'filled');
                    if showSamples
                        scatter3(points(:,1),points(:,2), points(:,3), [],values,'filled');
                    end
                    title({titlep1; ['Reconstructed Signal With ',num2str(numpnts),' Samples']});
                    colorbar();
                    view([az el])
                    xlabel('X1')
                    ylabel('X2')
                    zlabel('X3')
                    hold off
                end
            end
        case 4
            figure(1)
            clf(1)
            
            % Switch coordinate systems if necessary.
            [n1 n2 n3]=switcherfn(newpoints(:,1),newpoints(:,2),newpoints(:,3));
            newpoints(:,1) = n1;
            newpoints(:,2) = n2;
            newpoints(:,3) = n3;

            % Sorted list of times.
            times = sort(unique(newpoints(:,4)));
            
            % This varible determines which samples are shown in 
            % the current frame. Assuming the spacing between times 
            % is constant (which we know to be true since we 
            % generate them in generateSlicePoints.m) we can take 
            % the difference between the first two times and divide 
            % it by two.
            slack = (times(2)-times(1))/2;
            numframes = numel(times);
            
            % Find the bounds so the axes and color scale does not change at each
            % frame of the animation.
            xmin = min(newpoints(:,1));
            xmax = max(newpoints(:,1));
            ymin = min(newpoints(:,2));
            ymax = max(newpoints(:,2));
            zmin = min(newpoints(:,3));
            zmax = max(newpoints(:,3));
            bounds = [xmin xmax ymin ymax zmin zmax cbounds(1) cbounds(2)];
            
            if isVector
                % We have vector animations for this one case.
                if isAnimated
                    set(gca,'NextPlot','replaceChildren');
                    for j=1:numframes
                        title({titlep1; ['Time: ' num2str(times(j))]})
                        % Get all the new points for the current time.
                        ind=find(newpoints(:,4)==times(j));
                        easycone(newpoints(ind,1), newpoints(ind,2), newpoints(ind,3),...
                            guesses(ind,1),guesses(ind,2),guesses(ind,3),arrowsize);
                        view([az el])
                        pause(1/fps);
                    end
                else
                    plottingfn = @(tind) scatconehelper(newpoints,guesses,times,tind,arrowsize);
                    panelplotter(numframes,currp,plottingfn,times,az,el);
                end
            else
                if isAnimated
                    set(gca,'NextPlot','replaceChildren');
                    for j=1:numframes
                        % Get all the new points for the current time.
                        ind=find(newpoints(:,4)==times(j));
                        scatter3(newpoints(ind,1), newpoints(ind,2), newpoints(ind,3), 10, guesses(ind),'filled');
                        % Find any samples which are +/- slack distance from the current
                        % time in the animation and show them.
                        if showSamples
                            ind = find(points(:,4) > times(j)-slack & points(:,4) < times(j)+slack);
                            hold on
                            scatter3(points(ind,1),points(ind,2), points(ind,3), [],values(ind),'filled');
                            hold off
                        end
                        title({titlep1; ['Time: ' num2str(times(j))]})
                        axis(bounds);
                        view([az el])
                        pause(1/fps);
                    end
                else
                    plottingfn = @(tind) scatter3helper(newpoints,guesses,times,tind);
                    panelplotter(numframes,currp,plottingfn,times,az,el);
                end
            end
        otherwise
            disp('Dimension is too great to visualize');
    end
end

% This function is used throughout the main visualization code in the
% frame-by-frame option. numt is the number of times. currp is the current
% set of 4 frames.  plottingfn is used for the actual vis. times is a
% sorted list of all the times. az and el control the view angle.
function panelplotter(numt, currp, plottingfn, times,az,el)
% Number of 4 frame groups.
numpanels = ceil(numt/4);

% We're at the end of the frames.
if currp >= numpanels
    re = mod(numt,4);
    % If the number of frames is not divisible by four, then less than four
    % frames will be shown.
    switch re
        case 0
            subplot(2,2,1), plottingfn(numt-3); view([az el]); title(times(numt-3))
            subplot(2,2,2), plottingfn(numt-2); view([az el]); title(times(numt-2))
            subplot(2,2,3), plottingfn(numt-1); view([az el]); title(times(numt-1))
            subplot(2,2,4), plottingfn(numt); view([az el]); title(times(numt))
        case 1
            subplot(2,2,1), plottingfn(numt); view([az el]); title(times(numt))
        case 2
            subplot(2,2,1), plottingfn(numt-1); view([az el]); title(times(numt-1))
            subplot(2,2,2), plottingfn(numt); view([az el]); title(times(numt))
        case 3
            subplot(2,2,1), plottingfn(numt-2); view([az el]); title(times(numt-2))
            subplot(2,2,2), plottingfn(numt-1); view([az el]); title(times(numt-1))
            subplot(2,2,3), plottingfn(numt); view([az el]); title(times(numt))
    end
% We're at the beginning. Make sure there are 4 frames.
elseif currp <= 1
    if numt > 0
        subplot(2,2,1), plottingfn(1); view([az el]); title(times(1))
    end
    if numt > 1
        subplot(2,2,2), plottingfn(2); view([az el]); title(times(2))
    end
    if numt > 2
        subplot(2,2,3), plottingfn(3); view([az el]); title(times(3))
    end
    if numt > 3
        subplot(2,2,4), plottingfn(4); view([az el]); title(times(4))
    end
% We're somewhere in the middle.  We know there are 4 frames to show.
else
    subplot(2,2,1), plottingfn(4*currp-3); view([az el]); title(times(4*currp-3))
    subplot(2,2,2), plottingfn(4*currp-2); view([az el]); title(times(4*currp-2))
    subplot(2,2,3), plottingfn(4*currp-1); view([az el]); title(times(4*currp-1))
    subplot(2,2,4), plottingfn(4*currp); view([az el]); title(times(4*currp))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following functions are used in the panelplotter code for the       %
% frame-by-frame option.                                                  %
function easyconehelper(X1,X2,X3,g1,g2,g3,arrowsize,ind)                  %
easycone(X1(:,:,:,ind), X2(:,:,:,ind), X3(:,:,:,ind),g1(:,:,:,ind),...    %
    g2(:,:,:,ind), g3(:,:,:,ind), arrowsize);                             
                                                                          %
function easyslicehelper(X1,X2,X3,g,ind)
easyslice(X1(:,:,:,ind), X2(:,:,:,ind), X3(:,:,:,ind), g(:,:,:,ind));     %

function surfhelper(X1,X2,guesses,cbounds,ind)                            %
h = surf(X1(:,:,ind),X2(:,:,ind),guesses(:,:,ind));
set(h,'EdgeColor','none')                                                 %
title(['Time index: ' num2str(ind)])
caxis(cbounds)                                                            %

function plothelper(X1,guesses,ind)                                       %
plot(X1(:,ind),guesses(:,ind))
                                                                          %
function trisurfhelper(newpoints, guesses, times, timeind)
ind = find(newpoints(:,3)==times(timeind));
gridDelaunay = delaunay(newpoints(ind,1), newpoints(ind,2));              %
h = trisurf(gridDelaunay,newpoints(ind,1),newpoints(ind,2),guesses(ind));
set(h,'EdgeColor','none')                                                 %

function scatconehelper(newpoints, guesses, times, timeind, arrowsize)    %
ind = find(newpoints(:,4)==times(timeind));
easycone(newpoints(ind,1), newpoints(ind,2), newpoints(ind,3),...         %
    guesses(ind,1),guesses(ind,2),guesses(ind,3),arrowsize);
                                                                          %
function scatter3helper(newpoints, guesses, times, timeind)
ind = find(newpoints(:,4)==times(timeind));                               %
scatter3(newpoints(ind,1),newpoints(ind,2),newpoints(ind,3),10,...        %
    guesses(ind),'filled');                                               %
                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function switches the coordinate system of samples and new points.
function [X1p X2p X3p] = switchcoordsys(dispID, sampID, X1, X2, X3, pntdim)
if dispID == sampID
    X1p = X1;
    X2p = X2;
    X3p = X3;
    return;
end

switch sampID
    case 1 % Cartesian to...
        if dispID == 2 % polar.
            if pntdim == 2
                [X1p X2p] = cart2pol(X1,X2);
                X3p = 0;
            else
                [X1p X2p X3p] = cart2pol(X1,X2,X3);
            end
        elseif dispID == 3 % spherical.
            [X1p X2p X3p] = cart2sph(X1,X2,X3);
        end
    case 2 % Polar to...
        if dispID == 1 % cartesian.
            if pntdim == 2
                [X1p X2p] = pol2cart(X1,X2);
                X3p = 0;
            else
                [X1p X2p X3p] = pol2cart(X1,X2,X3);
            end
        elseif dispID == 3 % spherical.
            [X1p X2p X3p] = cart2sph(pol2cart(X1,X2,X3));
        end
    case 3 % Spherical to...
        if dispID == 1 % cartesian.
            [X1p X2p X3p] = sph2cart(X1,X2,X3);
        elseif dispID == 2 % polar
            [X1p X2p X3p] = cart2pol(sph2cart(X1,X2,X3));
        end
end

% function isvalid = checkbounds(bounds)
% i = 1;
% isvalid = 1;
% while i<numel(bounds)
%     if bounds(i) == bounds(i+1);
%         isvalid = 0;
%         return;
%     end
%     i = i+2;
% end
