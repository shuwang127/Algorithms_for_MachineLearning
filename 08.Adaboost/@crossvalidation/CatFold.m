function [Data, Labels] = CatFold(this, Data, Labels, N)

if(N > this.folds)
  error('N > total folds');
end

Data   = cat(2, Data, this.CrossDataSets{N}{1,1});
Labels = cat(2, Labels, this.CrossLabelsSets{N}{1,1});