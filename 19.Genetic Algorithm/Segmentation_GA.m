%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示遗传算法在计算机视觉中的应用
%实现如何利用遗传算法实现车牌图像的分割；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2015-04-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% read demo image
a=imread('testcard.jpg');  

% disp demo image
figure(1),imshow(a),title('Original Image')

% get dimension of image
[m,n]=size(a); 
a0=double(a);
h=1;                          
a1=zeros(m,n);

% Calculate 1D histogram 
for i=1:m
    for j=1:n
        for k=-h:h
            for w=-h:h;
                p=i+k;
                q=j+w;
                if (p<=0)|( p>m)
                    p=i;
                end
                if (q<=0)|(q>n)
                    q=j;
                end
                 a1(i,j)=a0(p,q)+a1(i,j);
             end
        end
        a2(i,j)=uint8(1/9*a1(i,j));
    end
end
fxy=zeros(256,256);

% Calculate 2D histogram 
for i=1:m
    for j=1:n
        c=a0(i,j);
        d=double(a2(i,j));
        fxy(c+1,d+1)=fxy(c+1,d+1)+1;
     end
end
Pxy=fxy/m/n;
P0=zeros(256,256);
Ui=zeros(256,256);
Uj=zeros(256,256);
P0(1,1)=Pxy(1,1);
for i=2:256
    P0(1,i)=P0(1,i-1)+Pxy(1,i);
end
for i=2:256
    P0(i,1)=P0(i-1,1)+Pxy(i,1);
end
for i=2:256
    for j=2:256
        P0(i,j)=P0(i,j-1)+P0(i-1,j)-P0(i-1,j-1)+Pxy(i,j);
    end
end
P1=ones(256,256)-P0;
Ui(1,1)=0;
for i=2:256
    Ui(1,i)=Ui(1,i-1)+(1-1)*Pxy(1,i);
end
for i=2:256
    Ui(i,1)=Ui(i-1,1)+(i-1)*Pxy(i,1);
end
for i=2:256
    for j=2:256
        Ui(i,j)=Ui(i,j-1)+Ui(i-1,j)-Ui(i-1,j-1)+(i-1)*Pxy(i,j);
    end
end
Uj(1,1)=0;
for i=2:256
    Uj(1,i)=Uj(1,i-1)+(i-1)*Pxy(1,i);
end
for i=2:256
    Uj(i,1)=Uj(i-1,1)+(1-1)*Pxy(i,1);
end
for i=2:256
    for j=2:256
        Uj(i,j)=Uj(i,j-1)+Uj(i-1,j)-Uj(i-1,j-1)+(j-1)*Pxy(i,j);
    end
end
uti=0;
utj=0;
for i=1:256
    for j=1:256
        uti=uti+(i-1)*Pxy(i,j); 
        utj=utj+(j-1)*Pxy(i,j);
    end
end

% Calculate max variance
hmax=0;
for i=1:256
    for j=1:256
       if P0(i,j)~=0&P1(i,j)~=0
        h(i,j)=((uti*P0(i,j)-Ui(i,j))^2+(utj*P0(i,j)-Uj(i,j))^2)/(P1(i,j)*P0(i,j));
       else
        h(i,j)=0;
       end
    end
end

hmax=max(h(:));
for i=1:256
    for j=1:256
        if h(i,j)==hmax
             s=i-1;
             t=j-1;
             continue;
         end
     end
end

% the result
z=ones(m,n);   
for i=1:m
    for j=1:n
      if a(i,j)<=s&a2(i,j)<=t 
       z(i,j)=0;      
      end
  end
end

% disp the result image
figure(2),imshow(z),title('GA-based car plate number segmentation');  