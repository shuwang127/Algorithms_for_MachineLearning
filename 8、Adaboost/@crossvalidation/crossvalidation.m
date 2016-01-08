function this = crossvalidation(folds)

if( folds == 1)
  error('folds should be >= 2');
end

this.folds = folds;

this.CrossDataSets   = cell(folds, 1);
this.CrossLabelsSets = cell(folds, 1);

this=class(this, 'crossvalidation') ;