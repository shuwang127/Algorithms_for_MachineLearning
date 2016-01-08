function code = TranslateToC (Learners, Weights, fid)

Weights = Weights ./ (sum(abs(Weights)));

fprintf(fid, ' %d\r\n ', length(Weights));

for i = 1 : length(Weights)
  Curr_Result = get_dim_and_tr(Learners{i});
  
  fprintf(fid, ' %f ', Weights(i));
  
  fprintf(fid, ' %d ', length(Curr_Result) / 3);
  
  for j = 1 : length(Curr_Result)
    fprintf(fid, ' %f ', Curr_Result(j));
  end
  
  fprintf(fid, '\r\n');
end

code = 1;