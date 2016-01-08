function [Data, Labels] = GetFold(this, N)

Data = [];
Labels = [];

if(N > this.folds)
  error('N > total folds');
end

Data   = cat(2, Data, this.CrossDataSets{N}{1,1});
Labels = cat(2, Labels, this.CrossLabelsSets{N}{1,1});
