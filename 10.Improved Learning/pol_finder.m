function policy=pol_finder(stat)

global NO_REPLICATIONS ITERMAX NA NS SMALL TPM TRM


for state=1:NS

[maxQfactor,index]=max(stat.Q(state,:));

policy(state)=index;

value_function(state)=maxQfactor;

end

policy

value_function

for state=1:NS

      for action=1:NA
      state
      action
      stat.Q(state,action)
      end

end





