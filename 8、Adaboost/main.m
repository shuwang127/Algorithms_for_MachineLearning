%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示Adaboost算法在计算机视觉中的应用
%基于Adaboost实现目标分类；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2013-09-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step1: reading Data from the file
file_data = load('Ionosphere.txt');
Data = file_data(:,1:end-1)';
Labels = file_data(:, end)';
Labels = Labels*2 - 1;

% boosting iterations
MaxIter = 100; 

% Step2: splitting data to training and control set
TrainData   = Data(:,1:2:end);
TrainLabels = Labels(1:2:end);

ControlData   = Data(:,2:2:end);
ControlLabels = Labels(2:2:end);

% Step3: constructing weak learner
weak_learner = tree_node_w(3); % pass the number of tree splits to the constructor

% Step4: training with Gentle AdaBoost
[RLearners RWeights] = RealAdaBoost(weak_learner, TrainData, TrainLabels, MaxIter);

% Step5: evaluating on control set
ResultR = sign(Classify(RLearners, RWeights, ControlData));

% Step6: calculating error
ErrorR  = sum(ControlLabels ~= ResultR) / length(ControlLabels)
