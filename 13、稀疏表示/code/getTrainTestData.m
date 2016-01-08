function [trainLabel, testLabel, trainSample, testSample,numClass] = ...
    getTrainTestData(datasetPath, fileExt, imgSide, percentTraining, normalizationFlag)

    imgSampleDir = dir([datasetPath, '/*.',fileExt]);
    
    numSample = length (imgSampleDir);
    sampleAll = zeros(2, numSample);
    
    imgResizedHeight = imgSide;
    imgResizedWidth  = imgSide;
    sampleImgAll = zeros(imgResizedHeight*imgResizedWidth, numSample);
    
    % Load all images

    for i = 1: numSample
        [sampleAll(1,i), sampleAll(2,i)] = getLabel(imgSampleDir(i).name);
        Img = imread([datasetPath, '/', imgSampleDir(i).name]);
        %load([datasetPath, '/', imgSampleDir(i).name]);
        if ndims(Img) == 3
        img = im2double(rgb2gray(Img));
        else
        img = im2double(Img);
        end
        img = adapthisteq(img);
        img = imresize(img, [imgResizedHeight, imgResizedWidth],'bicubic');
        sampleImgAll(:, i) = img(:);
    end
    
    % Normalization: normalize all samples with zero mean and unit variance
    if(normalizationFlag)
        sampleImgAllAvg=mean(sampleImgAll,1);
        sampleImgAll=sampleImgAll-ones(size(sampleImgAll,1),1)*sampleImgAllAvg;

%         for i = 1: size(sampleImgAll,2)
%             sampleImgAll(:,i) = sampleImgAll(:,i)/norm(sampleImgAll(:,i));
%         end
        sampleImgAllNorm = sqrt( sum(sampleImgAll.^2,1) );
        sampleImgAll = sampleImgAll./ (ones(size(sampleImgAll,1),1)*sampleImgAllNorm);
    end
        
    % Get index vectors for the training and testing samples    
    [trainIndex, testIndex] = getTrainTestInd(sampleAll, numSample, percentTraining);
    
    trainLabel = sampleAll(1,trainIndex);
    testLabel = sampleAll(1,testIndex);
    trainSample = sampleImgAll(:,trainIndex);
    testSample = sampleImgAll(:,testIndex);
        
    numClass = max(sampleAll(1,:));
   
end

function [label, num] =getLabel(filename)
    i = strfind(filename, '-');
    label = str2num(filename(1:i-1));
    j = strfind(filename(i+1:end),'.');
    num = str2num(filename(i+1:i+j));   
end

function [trainIndex, testIndex] = getTrainTestInd(sampleAll, numSample, percentTraining)
    
    trainIndex = zeros(1, numSample);
    
    numClass = max(sampleAll(2,:));
    for indClass = 1:numClass
        indSampleInClass=find(sampleAll(1,:)==indClass);
        numSampleInClass=length(indSampleInClass);
        numTrainSampleInClass= ceil(percentTraining*numSampleInClass);
        selectedIndexes = randperm(numSampleInClass);
        trainIndex(indSampleInClass(selectedIndexes(1:numTrainSampleInClass)))=1;
    end
    trainIndex = logical(trainIndex);
    testIndex = ~trainIndex;
end