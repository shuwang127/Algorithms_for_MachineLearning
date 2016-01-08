%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示增强学习算法在计算机视觉中的应用
%基于增强学习实现目标分类；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2014-02-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NO_REPLICATIONS ITERMAX NA NS SMALL TPM TRM LAMBDA

global3 % global parameters initialized 

stat=struct('Q',zeros(NS,NA),'iter',0,'old_action',1,'old_state',1,'current_state',1,'rimm',0,'total_reward',0);
done=0; % Pnemonic for simulation, 1 stands for end
        % 0 stands for continue 
while 0==done
    [stat,done]=jump_learn(stat);
end
   policy=pol_finder(stat);


