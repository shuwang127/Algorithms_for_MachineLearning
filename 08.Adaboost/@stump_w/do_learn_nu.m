%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.

function stump = do_learn_nu(stump, dataset, labels, weights)

%dataset=data_w(dataset) ;

%stump = set_distr(stump, get_sampl_weights(dataset)) ;

Distr = weights;

%[trainpat, traintarg] = get_train( dataset);
trainpat = dataset;
traintarg = labels;

tr_size = size(trainpat, 2);

T_MIN = zeros(3,size(trainpat,1));

for d = 1 : size(trainpat,1);

  [DS, IX] = sort(trainpat(d,:));

  TS = traintarg(IX);
  DiS = Distr(IX);
    
  lDS = length(DS);
  
  vPos = 0 * TS;
  vNeg = vPos;
  
  i = 1;
  j = 1;
  
  while i <= lDS
    k = 0;
    while i + k <= lDS && DS(i) == DS(i+k)
      if(TS(i+k) > 0)
        vPos(j) = vPos(j) + DiS(i+k);
      else
        vNeg(j) = vNeg(j) + DiS(i+k);
      end
      k = k + 1;
    end
    i = i + k;
    j = j + 1;
  end
  
  vNeg = vNeg(1:j-1);
  vPos = vPos(1:j-1);
  
  Error = zeros(1, j - 1);

  InvError = Error;
  
  IPos = vPos;
  INeg = vNeg;
  
  for i = 2 : length(IPos)
    IPos(i) = IPos(i-1) + vPos(i);
    INeg(i) = INeg(i-1) + vNeg(i);
  end
  
  Ntot = INeg(end);
  Ptot = IPos(end);
  
  for i = 1 : j - 1
    Error(i) = IPos(i) + Ntot - INeg(i);
    InvError(i) = INeg(i) + Ptot - IPos(i);
  end
  
  idx_of_err_min = find(Error == min(Error));
  if(length(idx_of_err_min) < 1)
      idx_of_err_min = 1;  
  end
  
  if(length(idx_of_err_min) <1)
    idx_of_err_min = idx_of_err_min;
  end
  idx_of_err_min = idx_of_err_min(1);
  
  idx_of_inv_err_min = find(InvError == min(InvError));
  
  if(length(idx_of_inv_err_min) < 1)
      idx_of_inv_err_min = 1;  
  end
  
  idx_of_inv_err_min = idx_of_inv_err_min(1);
  
  if(Error(idx_of_err_min) < InvError(idx_of_inv_err_min))
    T_MIN(1,d) = Error(idx_of_err_min);
    T_MIN(2,d) = idx_of_err_min;
    T_MIN(3,d) = -1;
  else
    T_MIN(1,d) = InvError(idx_of_inv_err_min);
    T_MIN(2,d) = idx_of_inv_err_min;
    T_MIN(3,d) = 1;
  end
  
end

best_dim = find(T_MIN(1,:) == min(T_MIN(1,:)));
stump.t_dim = best_dim(1);

TDS = sort(trainpat(stump.t_dim,:));

lDS = length(TDS);

DS = TDS * 0;

i = 1;
j = 1;

while i <= lDS
  k = 0;
  while i + k <= lDS && TDS(i) == TDS(i+k) 
    DS(j) = TDS(i);
    k = k + 1;
  end
  i = i + k;
  j = j + 1;
end

DS = DS(1:j-1);

stump.threshold = (DS(T_MIN(2,stump.t_dim)) + DS(min(T_MIN(2,stump.t_dim) + 1, length(DS)))) / 2;
stump.signum = T_MIN(3,stump.t_dim);




