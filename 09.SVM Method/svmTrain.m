function svm = svmTrain(svmType,X,Y,ker,p1,p2)
% SVM Classification:
%   svm = svmTrain('svc_c',x,y,ker,C); 
%   svm = svmTrain('svc_nu',x,y,ker,nu); 
%
% One-Class SVM:
%   svm = svmTrain('svm_one_class',x,[],ker,nu);
%
% SVM Regression:
%   svm = svmTrain('svr_epsilon',x,y,ker,C,e); 
%   svm = svmTrain('svr_nu',x,y,ker,C,nu); 

options = optimset;
options.LargeScale = 'off';
options.Display = 'off';

switch svmType
    case 'svc_c',
        
        C = p1;

        n = length(Y);
        H = (Y'*Y).*kernel(ker,X,X);

        f = -ones(n,1);
        A = [];
        b = [];
        Aeq = Y;
        beq = 0;
        lb = zeros(n,1);
        ub = C*ones(n,1);
        a0 = zeros(n,1);
        
        [a,fval,eXitflag,output,lambda]  = quadprog(H,f,A,b,Aeq,beq,lb,ub,a0,options);                    
    
    case 'svc_nu',
        
        nu = p1;
       
        n = length(Y);
        H = (Y'*Y).*kernel(ker,X,X);

        f = zeros(n,1);
        A = -ones(1,n);
        b = -nu;
        Aeq = Y;
        beq = 0;
        lb = zeros(n,1);
        ub = ones(n,1)/n;
        a0 = zeros(n,1);
        
        [a,fval,eXitflag,output,lambda]  = quadprog(H,f,A,b,Aeq,beq,lb,ub,a0,options);                    

    case 'svm_one_class',
        
        nu = p1;
    
        n = size(X,2);
        H = kernel(ker,X,X);

        f = zeros(n,1);
        for i = 1:n
            f(i,:) = -kernel(ker,X(:,i),X(:,i));
        end
        A = [];
        b = [];
        Aeq = ones(1,n);
        beq = 1;
        lb = zeros(n,1);
        ub = ones(n,1)/(nu*n);
        a0 = zeros(n,1);

        [a,fval,eXitflag,output,lambda]  = quadprog(H,f,A,b,Aeq,beq,lb,ub,a0,options);                    
        
    case 'svr_epsilon',
        
        C = p1;
        e = p2;
        
        n = length(Y);
        Q = kernel(ker,X,X);
        H = [Q,-Q;-Q,Q];
        f = [e*ones(n,1)-Y';e*ones(n,1)+Y'];         
        %f = [e*ones(n,1)+Y';e*ones(n,1)-Y'];
        A = [];
        b = [];
        Aeq = [ones(1,n),-ones(1,n)];
        beq = 0;
        lb = zeros(2*n,1);               
        ub = C*ones(2*n,1);
        a0 = zeros(2*n,1);
        
        [a,fval,eXitflag,output,lambda]  = quadprog(H,f,A,b,Aeq,beq,lb,ub,a0,options);  
        a = a(1:n)-a(n+1:end);

    case 'svr_nu',

        C = p1;
        nu = p2;

        n = length(Y);
        Q = kernel(ker,X,X);
        H = [Q,-Q;-Q,Q];
        f = [-Y';+Y'];         
        %f = [+Y';-Y'];
        A = [];
        b = [];
        Aeq = [ones(1,n),-ones(1,n);ones(1,2*n)];
        beq = [0;C*n*nu];
        lb = zeros(2*n,1);               
        ub = C*ones(2*n,1);
        a0 = zeros(2*n,1);
        
        [a,fval,eXitflag,output,lambda]  = quadprog(H,f,A,b,Aeq,beq,lb,ub,a0,options);            
        a = a(1:n)-a(n+1:end);

    otherwise,
end


eXitflag

% output svm

svm.type = svmType;
svm.ker = ker;
svm.x = X;
svm.y = Y;
svm.a = a';