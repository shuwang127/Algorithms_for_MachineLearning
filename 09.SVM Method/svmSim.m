function Yd = svmSim(svm,Xt)

% parameters
cathe = 10e+6;                
nx = size(svm.x,2);           
nt = size(Xt,2);              
block = ceil(nx*nt/cathe);    
num = ceil(nt/block);          

for i = 1:block
    if (i==block)
        index = [(i-1)*num+1:nt];
    else
        index = (i-1)*num+[1:num];
    end

    Yd(index) = svmSim_block(svm,Xt(:,index));         
end

% ------------------------------------------------------------%

function Yd = svmSim_block(svm,Xt);

type = svm.type;
ker = svm.ker;
X = svm.x;
Y = svm.y;
a = svm.a;

% test output
epsilon = 1e-8;                 
i_sv = find(abs(a)>epsilon);         

switch type
    case 'svc_c',        
        tmp = (a.*Y)*kernel(ker,X,X(:,i_sv));         
        b = Y(i_sv)-tmp;
        b = mean(b);
        tmp =  (a.*Y)*kernel(ker,X,Xt);
        tmp = tmp+b;
        Yd = sign(tmp);
        
    case 'svc_nu', 
        tmp = (a.*Y)*kernel(ker,X,X(:,i_sv));          
        b = Y(i_sv)-tmp;
        b = mean(b);
        tmp =  (a.*Y)*kernel(ker,X,Xt);
        Yd = sign(tmp+b);
        
    case 'svm_one_class',        
        n_sv = length(i_sv);
        tmp1 = zeros(n_sv,1);
        for i = 1:n_sv
            index = i_sv(i);
            tmp1(i) = kernel(ker,X(:,index),X(:,index));
        end

        tmp2 = 2*a*kernel(ker,X,X(:,i_sv));           
        tmp3 = sum(sum(a'*a.*kernel(ker,X,X)));    

        R_square = tmp1-tmp2'+tmp3;
        R_square = mean(R_square);                       

        nt = size(Xt,2);                

        tmp4 = zeros(nt,1);             
        for i = 1:nt
            tmp4(i) = kernel(ker,Xt(:,i),Xt(:,i));
        end
    
        tmp5 = 2*a*kernel(ker,X,Xt);               
        Yd = sign(tmp4-tmp5'+tmp3-R_square);

    case 'svr_epsilon',
        
        tmp = a*kernel(ker,X,X(:,i_sv));  
        b = Y(i_sv)-tmp;                  
        %b = Y(i_sv)+tmp;
        b = mean(b);

        tmp =  a*kernel(ker,X,Xt);         
        %tmp =  -a*kernel(ker,X,Xt);
        Yd = (tmp+b);        
        
    case 'svr_nu',
      
        tmp = a*kernel(ker,X,X(:,i_sv));   
        b = Y(i_sv)-tmp;                   
        %b = Y(i_sv)+tmp;
        b = mean(b);

        tmp =  a*kernel(ker,X,Xt);         
        %tmp =  -a*kernel(ker,X,Xt);
        Yd = (tmp+b);        
        
    otherwise,
end