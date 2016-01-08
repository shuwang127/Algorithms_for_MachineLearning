function this = Initialize(this, Data, Labels)

if( size(Labels,1) ~= 1 || length(size(Labels)) ~= 2)
  error('Labels should be a (1,N) matrix.');
end

if( size(Labels,2) ~= size(Data,2))
  error('Data should be (M,N) matrix and Labels (1,N)');
end

for i = 1 : length(Data)
    fold_i = floor(rand * (this.folds - eps) + 1);
    if( i > this.folds )
        this.CrossDataSets{fold_i} = {cat(2, Data(:, i), this.CrossDataSets{fold_i}{1,1})};
        this.CrossLabelsSets{fold_i} = {cat(2, Labels(:, i), this.CrossLabelsSets{fold_i}{1,1})};
    else
        this.CrossDataSets{i} = {cat(2, Data(:, i), this.CrossDataSets{i})};
        this.CrossLabelsSets{i} = {cat(2, Labels(:, i), this.CrossLabelsSets{i})};
    end
end
