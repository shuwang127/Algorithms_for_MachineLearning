function stat=qlearn(stat)

global NO_REPLICATIONS ITERMAX NA NS SMALL TPM TRM LAMBDA

% Q-Learning 

% Finding the Max factor in the current state 

q_next=max(stat.Q(stat.current_state,:));


stat.iter=stat.iter+1;


%learn_rate=1/(stat.iter);

learn_rate=log(stat.iter+1)/(stat.iter+1);

%learn_rate=0.5*300/(300+stat.iter);

q=stat.Q(stat.old_state,stat.old_action);

q=q*(1-learn_rate)+(learn_rate*(stat.rimm+(LAMBDA*q_next)));

stat.Q(stat.old_state,stat.old_action)=q;




