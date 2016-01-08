file_data = load('Diabetes.txt');
Data = file_data(:,1:end-1)';
Labels = file_data(:, end)';
Labels = Labels*2 - 1;

% Data = Data';
% Labels = Labels';

weak_learner = tree_node_w(2);

% Step1: training with Gentle AdaBoost
[RLearners RWeights] = RealAdaBoost(weak_learner, Data, Labels, 200);

% Step2: training with Modest AdaBoost
[MLearners MWeights] = ModestAdaBoost(weak_learner, Data, Labels, 200);

fid = fopen('RealBoost.txt','w');
TranslateToC(RLearners, RWeights, fid);
fclose(fid);

fid = fopen('ModestBoost.txt','w');
TranslateToC(MLearners, MWeights, fid);
fclose(fid);