function action=action_selector(stat)

global NA 


ran=rand(1);

candidate=1;

sum=1/NA;

complete=0;

% Selecting each action with equal probability 

while 0==complete

        if ran<sum
        % action selected 
        action=candidate;
        complete=1;
        
        else          
        % test if ran is associated with next action 
        candidate=candidate+1;
        sum=sum+(1/NA);
        end
end

