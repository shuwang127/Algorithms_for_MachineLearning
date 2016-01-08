%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示EM算法在计算机视觉中的应用
%基于EM算法实现目标分类；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2013-09-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1=[0 6];
mu2=[0 0];
mu3=[6,-6];
sigma1=[8 0; 0 1]
r1=mvnrnd(mu1,sigma1,100);
r2=mvnrnd(mu2,sigma1,100);
r3=mvnrnd(mu3,sigma1,100)
data(1:100,:)=r1;
data(101:200,:)=r2;
data(201:300,:)=r3;
[n,d]=size(data)
M=3;
sigma=zeros(d,d,M);
Mu=zeros(M,d);
size(Mu)
priors=zeros(1,M);
[priors, Mu,sigma] = EM_init_kmeans(data', M)

max_iter=20;
for i=1:max_iter
    p=zeros(M,n);%%p=p(x|j)
    E=zeros(M,n);
    %% E-step
    for j=1:M
       if det(sigma(:,:,j))==0
           sigma(:,:,j)=ones*exp(-100);%保证矩阵可逆
           
       end
       detS=det(sigma(:,:,j));
       invS=inv(sigma(:,:,j));
       for k=1:n         
           p(j,k)=(2*pi)^(-d/2)*(detS)^(-1/2)*exp(-1/2*(data(k,:)-Mu(j,:))*invS*(data(k,:)-Mu(j,:))');
           E(j,k)=p(j,k)*priors(j);        
       end
    end
    E=E./(ones(M,1)*sum(E));
    %% M-step
    w=zeros(1,M);
    Mu2=zeros(M,d);
    sigma2=zeros(d,d,M);
    % update mean 
    for k=1:M
        for j=1:n
            w(k)=w(k)+E(k,j);
            Mu2(k,:)=Mu2(k,:)+E(k,j)*data(j,:);       
        end
        Mu2(k,:)=Mu2(k,:)./w(k);
    end
    % update covariance
    for k=1:M
        for j=1:n
            dist=(Mu2(k,:)-data(j,:))*(Mu2(k,:)-data(j,:))';
            sigma2(:,:,k)=sigma2(:,:,k)+E(k,j)*dist;            
        end
        sigma2(:,:,k)=sigma2(:,:,k)./w(k);
        sigma2(:,:,k)=diag(diag(sigma2(:,:,k)));
        
    end
    %update prior probability
    priors=w./n;
    Mu=Mu2;
    sigma=sigma2;
    
end
%% data classification
[x,y]=meshgrid(-10:.1:10,-10:.1:10);
figure(2)
subplot(1,2,1)
plot(r1(:,1),r1(:,2),'*')
hold on
plot(r2(:,1),r2(:,2),'o')
hold on
plot(r3(:,1),r3(:,2),'+')
title('原始数据');
xlim([-10 10])
ylim([-10 10])

z1=mvnpdf([x(:) y(:)],mu1,sigma1);
z=zeros(size(z1));
z1 = reshape(z1,length(x),length(y));
hold on
contour(x,y,z1,[0.11 0.1 0.08 0.06 0.04 0.02 0.01 0.001 0.0001]);
z2=mvnpdf([x(:) y(:)],mu2,sigma1);
z2 = reshape(z2,length(x),length(y));
hold on
contour(x,y,z2,[0.11 0.1 0.08 0.06 0.04 0.02 0.01 0.001 0.0001]);
z3=mvnpdf([x(:) y(:)],mu3,sigma1);
z3 = reshape(z3,length(x),length(y));
hold on
contour(x,y,z3,[0.11 0.1 0.08 0.06 0.04 0.02 0.01 0.001 0.0001]);
subplot(1,2,2)

C=['*','O','+'];

maxE=max(E);
result=zeros(n,d,M);
%z=zeros(size(x,1)*size(y,1),1);
for i=1:M
    in=(E(i,:)==maxE);
    acount=find(E(i,:)==maxE);
    number=size(acount,2)
    if ~isempty(acount)
    result=zeros(size(acount,2),d);
    m=1
    x=acount(m);
    %m=1;
    for j=1:n
        if j==x;
            
            result(m,:)=data(x,:);
            m=m+1;
            if m>size(acount,2)
              break;  
            end
            x=acount(m);
            
            
        end
    end
  plot(result(:,1),result(:,2),C(i));
    end
  %x=diag(sigma(:,:,i))
  [x,y]=meshgrid(-10:.1:10,-10:.1:10);

   z=z+mvnpdf([x(:) y(:)],Mu(i,:),sigma(:,:,i))/priors(i);
  hold on
  
end
z = reshape(z,length(x),length(y));
contour(x,y,z,[0.025 0.02 0.01]);
title('分类后数据')

 
