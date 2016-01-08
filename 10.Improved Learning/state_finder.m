function candidate=state_finder(stat)

global NO_REPLICATIONS ITERMAX NA NS SMALL TPM TRM

ran=rand(1);

old_action=stat.old_action;
old_state=stat.old_state;

sum=TPM(old_state,1,old_action);

candidate=1;

complete=0;

while 0==complete

        if ran<sum 
        
        complete=1;
        
        else
        
        candidate=candidate+1;

        sum=sum+TPM(old_state,candidate,old_action);

        end
        
end
                  
