%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示Kmeans聚类算法在计算机视觉中的应用
%实现如何利用Kmeans聚类实现图像的分割；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2014-10-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kmeans_demo1()
clear;close all;clc;
%% 读取测试图像
im = imread('city.jpg');
imshow(im), title('Imput image');
%% 转换图像的颜色空间得到样本
cform = makecform('srgb2lab');
lab = applycform(im,cform);
ab = double(lab(:,:,2:3));
nrows = size(lab,1); ncols = size(lab,2);
X = reshape(ab,nrows*ncols,2)';
figure, scatter(X(1,:)',X(2,:)',3,'filled');  box on; %显示颜色空间转换后的二维样本空间分布
%print -dpdf 2D1.pdf
%% 对样本空间进行Kmeans聚类
k = 5; % 聚类个数
max_iter = 100; %最大迭代次数

[centroids, labels] = run_kmeans(X, k, max_iter); 

%% 显示聚类分割结果
figure, scatter(X(1,:)',X(2,:)',3,labels,'filled'); %显示二维样本空间聚类效果
hold on; scatter(centroids(1,:),centroids(2,:),60,'r','filled')
hold on; scatter(centroids(1,:),centroids(2,:),30,'g','filled')
box on; hold off;
%print -dpdf 2D2.pdf

pixel_labels = reshape(labels,nrows,ncols);
rgb_labels = label2rgb(pixel_labels);
figure, imshow(rgb_labels), title('Segmented Image');
%print -dpdf Seg.pdf
end

function [centroids, labels] = run_kmeans(X, k, max_iter)
% 该函数实现Kmeans聚类
% 输入参数：
%                   X为输入样本集，dxN
%                   k为聚类中心个数
%                   max_iter为kemans聚类的最大迭代的次数
% 输出参数：
%                   centroids为聚类中心 dxk
%                   labels为样本的类别标记

%% 采用K-means++算法初始化聚类中心
  centroids = X(:,1+round(rand*(size(X,2)-1)));
  labels = ones(1,size(X,2));
  for i = 2:k
        D = X-centroids(:,labels);
        D = cumsum(sqrt(dot(D,D,1)));
        if D(end) == 0, centroids(:,i:k) = X(:,ones(1,k-i+1)); return; end
        centroids(:,i) = X(:,find(rand < D/D(end),1));
        [~,labels] = max(bsxfun(@minus,2*real(centroids'*X),dot(centroids,centroids,1).'));
  end
  
%% 标准Kmeans算法
  for iter = 1:max_iter
        for i = 1:k, l = labels==i; centroids(:,i) = sum(X(:,l),2)/sum(l); end
        [~,labels] = max(bsxfun(@minus,2*real(centroids'*X),dot(centroids,centroids,1).'),[],1);
  end
  
end


