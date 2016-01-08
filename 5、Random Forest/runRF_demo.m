path = ['E:\works\book\7£¨»úÆ÷Ñ§Ï°20½²£©\Code\5¡¢Random Forest\'];
trainData = [path 'satimage.tra'];
testData = [path 'satimage.txt'];
classwt=0;
cat0 = 0;
runParam = [6 1 50 10 1 0];
impOpt = [0 0 0];
proCom = [0 4435 0 0 0];
missingVal=0;
saveForest=[1 0 0];
runForest=[0 0];
outParam = [1,0,0,0,0,0,0,0,0,0];
msm = 0;
seed = 4351;

[errtr, errts, prox, trees, predictts, varimp, scale] = ...
  runRF(trainData, testData, classwt, cat0, msm, runParam, ...
    impOpt, proCom, missingVal, saveForest, runForest, outParam, seed);