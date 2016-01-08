function main

clc;clear all;close all;
img=imread('testcard.jpg');
defaultConfig.xy          = 10*rand(50,2);
defaultConfig.dmat        = [];
defaultConfig.popSize     = 100;
defaultConfig.numIter     = 1e4;
defaultConfig.showProg    = true;
defaultConfig.showResult  = true;
defaultConfig.showWaitbar = false;
% define the position of the plate number;
aa=[157,226];bb=[140,140];cc=[163,163];dd=[141,163];ee=[157,157];ff=[226,226];
% Interpret user configuration inputs
    if ~nargin
        userConfig = struct();
    elseif isstruct(varargin{1})
        userConfig = varargin{1};
    else
        try
            userConfig = struct(varargin{:});
        catch
            error('Expected inputs are either a structure or parameter/value pairs');
        end
    end
    
    % Override default configuration with user inputs
    configStruct = get_config(defaultConfig,userConfig);
    
    % Extract configuration
    xy          = configStruct.xy;
    dmat        = configStruct.dmat;
    popSize     = configStruct.popSize;
    numIter     = configStruct.numIter;
    showProg    = configStruct.showProg;
    showResult  = configStruct.showResult;
    showWaitbar = configStruct.showWaitbar;
    if isempty(dmat)
        nPoints = size(xy,1);
        a = meshgrid(1:nPoints);
        dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nPoints,nPoints);
    end
    
    % Verify Inputs
    [N,dims] = size(xy);
    [nr,nc] = size(dmat);
    if N ~= nr || N ~= nc
        error('Invalid XY or DMAT inputs!')
    end
    n = N;
    
    % Sanity Checks
    popSize     = 4*ceil(popSize/4);
    numIter     = max(1,round(real(numIter(1))));
    showProg    = logical(showProg(1));
    showResult  = logical(showResult(1));
    showWaitbar = logical(showWaitbar(1));
    
    % Initialize the Population
    pop = zeros(popSize,n);
    pop(1,:) = (1:n);
    for k = 2:popSize
        pop(k,:) = randperm(n);
    end
    
    % Run the GA
    globalMin = Inf;
    totalDist = zeros(1,popSize);
    distHistory = zeros(1,numIter);
    tmpPop = zeros(4,n);
    newPop = zeros(popSize,n);
    if showProg
%         figure('Name','TSP_GA | Current Best Solution','Numbertitle','off');
        hAx = gca;
    end
    if showWaitbar
        hWait = waitbar(0,'Searching for near-optimal solution ...');
    end
    for iter = 1:numIter
        % Evaluate Each Population Member (Calculate Total Distance)
        for p = 1:popSize
            d = dmat(pop(p,n),pop(p,1)); % Closed Path
            for k = 2:n
                d = d + dmat(pop(p,k-1),pop(p,k));
            end
            totalDist(p) = d;
        end
        
        % Find the Best Route in the Population
        [minDist,index] = min(totalDist);
        distHistory(iter) = minDist;
        if minDist < globalMin
            globalMin = minDist;
            optRoute = pop(index,:);
            if showProg
                % Plot the Best Route
                rte = optRoute([1:n 1]);
%                 if dims > 2, plot3(hAx,xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
%                 else plot(hAx,xy(rte,1),xy(rte,2),'r.-'); end
%                 title(hAx,sprintf('Total Distance = %1.4f, Iteration = %d',minDist,iter));
                drawnow;
            end
        end
        
        % Genetic Algorithm Operators
        randomOrder = randperm(popSize);
        for p = 4:4:popSize
            rtes = pop(randomOrder(p-3:p),:);
            dists = totalDist(randomOrder(p-3:p));
            [ignore,idx] = min(dists); %#ok
            bestOf4Route = rtes(idx,:);
            routeInsertionPoints = sort(ceil(n*rand(1,2)));
            I = routeInsertionPoints(1);
            J = routeInsertionPoints(2);
            for k = 1:4 % Mutate the Best to get Three New Routes
                tmpPop(k,:) = bestOf4Route;
                switch k
                    case 2 % Flip
                        tmpPop(k,I:J) = tmpPop(k,J:-1:I);
                    case 3 % Swap
                        tmpPop(k,[I J]) = tmpPop(k,[J I]);
                    case 4 % Slide
                        tmpPop(k,I:J) = tmpPop(k,[I+1:J I]);
                    otherwise % Do Nothing
                end
            end
            newPop(p-3:p,:) = tmpPop;
        end
        pop = newPop;
        
        % Update the waitbar
        if showWaitbar && ~mod(iter,ceil(numIter/325))
            waitbar(iter/numIter,hWait);
        end
        
    end
    if showWaitbar
        close(hWait);
    end
%     if showResult
%         % Plots the GA Results
% %         figure('Name','TSP_GA | Results','Numbertitle','off');
% %         subplot(2,2,1);
%         pclr = ~get(0,'DefaultAxesColor');
%         if dims > 2, plot3(xy(:,1),xy(:,2),xy(:,3),'.','Color',pclr);
%         else plot(xy(:,1),xy(:,2),'.','Color',pclr); end
%         title('City Locations');
%         subplot(2,2,2);
%         imagesc(dmat(optRoute,optRoute));
%         title('Distance Matrix');
%         subplot(2,2,3);
%         rte = optRoute([1:n 1]);
%         if dims > 2, plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
%         else plot(xy(rte,1),xy(rte,2),'r.-'); end
%         title(sprintf('Total Distance = %1.4f',minDist));
%         subplot(2,2,4);
%         plot(distHistory,'b','LineWidth',2);
%         title('Best Solution History');
%         set(gca,'XLim',[0 numIter+1],'YLim',[0 1.1*max([1 distHistory])]);
%     end
    imshow(img);hold on;
    plot(aa,bb,'w.-',aa,cc,'w.-',ee,dd,'w.-',ff,dd,'w.-');
    % Return Output
    if nargout
        resultStruct = struct( ...
            'xy',          xy, ...
            'dmat',        dmat, ...
            'popSize',     popSize, ...
            'numIter',     numIter, ...
            'showProg',    showProg, ...
            'showResult',  showResult, ...
            'showWaitbar', showWaitbar, ...
            'optRoute',    optRoute, ...
            'minDist',     minDist);
        
        varargout = {resultStruct};
    end
    figure(2);
   subimg=img(dd(1):dd(2),aa(1):aa(2)); imshow(subimg);
end
% Subfunction to override the default configuration with user inputs
function config = get_config(defaultConfig,userConfig)
    
    % Initialize the configuration structure as the default
    config = defaultConfig;
    
    % Extract the field names of the default configuration structure
    defaultFields = fieldnames(defaultConfig);
    
    % Extract the field names of the user configuration structure
    userFields = fieldnames(userConfig);
    nUserFields = length(userFields);
    
    % Override any default configuration fields with user values
    for i = 1:nUserFields
        userField = userFields{i};
        isField = strcmpi(defaultFields,userField);
        if nnz(isField) == 1
            thisField = defaultFields{isField};
            config.(thisField) = userConfig.(userField);
        end
    end
    
end