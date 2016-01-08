function [result2, relative_error2] = src(trn,tst, solver)
repeat = 1;
for r = 1:repeat
% Normalising the training set
for i=1:size(trn.X,2)
      Atilde(:,i)=trn.X(:,i)/norm(trn.X(:,i),2);
end

n = size(trn.X,2);    % total number of samples in the training database
k = max(trn.y);       % number of classes
onesvecs = full(ind2vec(trn.y));    % preparing classification target
n_test = size(tst.X,2);    % number of test data to be classified
if strcmp(solver, 'omp')
   param.L = 9;   % parameter for Sparse Solver
   param.eps=0.1;
end
if strcmp(solver, 'apg')
   param.regul='l1';   % parameter for Sparse Solver
   param.loss='square';
   beta0=zeros(size(Atilde,2),1);
   param.lambda = 0.07;
end
% start classifcation
for iter = 1:n_test    % looping through all the test samples
      ytilde = tst.X(:,iter);
      ytilde = ytilde/norm(ytilde);        % normalizing the test sample
      if strcmp(solver, 'omp'); xp = mexOMP(ytilde,Atilde,param);end
      if strcmp(solver, 'apg');  xp =mexFistaFlat(ytilde,Atilde,beta0,param);end
      % decision making
      for i=1:k
            deltavec(:,i) = onesvecs(i,:)'.* xp;
            residual2(r,iter,i) = norm(ytilde-Atilde*deltavec(:,i));
      end
      [B,IX] = sort(residual2);
      result(iter) = IX(1);
end
end

res2 = sum(residual2,1);

for iter = 1:n_test
    [B, Ind] = sort(res2(1,iter,:));
    result2(iter) = Ind(1);
end

% finding the relative classification error
relative_error2 = size(find(abs(tst.y-result2)),2)/size(tst.y,2);