function [errtr, errts, prox, trees, predictts, varimp, scale] = ...
  runRF(trainData, testData, classwt, cat0, msm, runParam, ...
  impOpt, proCom, missingVal, saveForest, runForest, outParam, seed)
  
  if (trainData == testData) 
    data = textread(trainData);
    [m, n] = size(data);
    randomNum = randperm(m);
    randomData = data(randomNum, :);
    ntest = ceil(m / 3);
    ntrain = m - ntest;
    x = (randomData(1 : ntrain, 1 : n - 1));
    x = x';
    cl = (randomData(1 : ntrain, n));
    xts = (data(ntrain + 1 : m, 1 : n - 1));
    xts = xts';
    clts = (data(ntrain + 1 : m, n));
  else
    data1 = textread(trainData);
    data2 = textread(testData);
    [m1, n1] = size(data1);
    [m2, n2] = size(data2); 
    ntest = m2;
    ntrain = m1;
    x = (data1(1 : ntrain, 1 : n1 - 1));
    x = x';
    cl = (data1(1 : ntrain, n1));
    xts = (data2(1 : ntest, 1 : n2 - 1));
    xts = xts';
    clts = (data2(1 : ntest, n2));
  end

  if (nargin < 13) 
    seed = 4351;
  end
  if (nargin < 12) 
    outParam = 0;
  end
  if (nargin < 11) 
    runForest = 0;
  end
  if (nargin < 10) 
    saveForest = 0;
  end
  if (nargin < 9) 
    missingVal = 0;
  end
  if (nargin < 8) 
    proCom = 0;
  end
  if (nargin < 7) 
    impOpt = 0;
  end
  if (nargin < 6) 
    runParam = 0;
  end
  if (nargin < 5) 
    msm = 0;
  end
  if (nargin < 4) 
    cat0 = 0;
  end
  if (nargin < 3) 
    classwt = 0;
  end
  if (nargin < 2) 
    xts = 0;
    clts = 0;
  end
  uniquecl = unique(cl);
  nclass = size(uniquecl, 1) * size(uniquecl, 2);
  remain = trainData;
  while true
    [dataName, remain] = strtok(remain,'\');
    if isempty(remain)
      break;
    end
  end
  dataName = strtok(dataName, '.');

  x = single(x);        %get train x
  cl = int32(cl);           %get train label
  xts = single(xts);      %get test x
  clts = int32(clts);     %get test label
  classwt = single(classwt);
  cat0 = int32(cat0);
  msm = int32(msm);
  runParam = int32(runParam);
  impOpt = int32(impOpt);
  proCom = int32(proCom);
  missingVal = single(missingVal);
  saveForest = int32(saveForest);
  runForest = int32(runForest);
  outParam = int32(outParam);
  seed = int32(seed);

  [errtr, errts, prox, trees, predictts, varimp, scale] = ...
    RF(nclass, x, cl, xts, clts, classwt, cat0, msm, runParam, impOpt, proCom, ...
      missingVal, saveForest, runForest, outParam, seed, dataName);
  errtr = errtr';
  errts = errts';
  prox = prox';
  predictts = predictts';
  varimp = varimp';
  scale = scale';
end