function Result = Classify(Learners, Weights, Data)

Result = zeros(1, size(Data, 2));

for i = 1 : length(Weights)
  lrn_out = calc_output(Learners{i}, Data);
  Result = Result + lrn_out * Weights(i);
end