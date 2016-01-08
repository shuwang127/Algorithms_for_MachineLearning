%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示KNN算法在计算机视觉中的应用
%实现如何利用KNN算法进行聚类；
%Modi: NUDT-VAP
%环境：Win7，Matlab2012b
%时间：2014-10-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main
trainData = [
    0.6213    0.5226    0.9797    0.9568    0.8801    0.8757    0.1730    0.2714    0.2523
    0.7373    0.8939    0.6614    0.0118    0.1991    0.0648    0.2987    0.2844    0.4692
    ];
trainClass = [
    1     1     1     2     2     2     3     3     3
    ];
testData = [
    0.9883    0.5828    0.4235    0.5155    0.3340
    0.4329    0.2259    0.5798    0.7604    0.5298
    ];

% main
testClass = cvKnn(testData, trainData, trainClass);

% plot prototype vectors
classLabel = unique(trainClass);
nClass     = length(classLabel);
plotLabel = {'r*', 'g*', 'b*'};
figure;
for i=1:nClass
    A = trainData(:, trainClass == classLabel(i));
    plot(A(1,:), A(2,:), plotLabel{i});
    hold on;
end

% plot classifiee vectors
plotLabel = {'ro', 'go', 'bo'};
for i=1:nClass
    A = testData(:, testClass == classLabel(i));
    plot(A(1,:), A(2,:), plotLabel{i});
    hold on;
end
legend('1: prototype','2: prototype', '3: prototype', '1: classifiee', '2: classifiee', '3: classifiee', 'Location', 'NorthWest');
title('K nearest neighbor');
hold off;