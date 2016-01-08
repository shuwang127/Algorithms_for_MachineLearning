function [stat,done]=jump_learn(stat)

global NO_REPLICATIONS ITERMAX NA NS SMALL TPM TRM

% This function simulates a jump and also updates the learning stats

old_state=stat.old_state;

old_action=stat.old_action;

% Determine current state 

current_state=state_finder(stat);

% Record Feedback in stat 

stat.current_state=current_state;

stat.rimm=TRM(old_state,current_state,old_action);

% DO LEARNING 

stat=qlearn(stat);

% Select next action 

next_action=action_selector(stat);

% Get ready to get out of this function 

stat.old_state=current_state;

stat.old_action=next_action;

	if stat.iter>=ITERMAX
	
	% Learning should end 

	done=1;
	
	else
	
	done=0;
      
        end







