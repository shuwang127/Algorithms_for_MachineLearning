function [Rdx,RdX,RdXt,Rdy,RdY,RdYt] = plsrd(E0,F0,T,h)
%E0-标准化后的自变量数据
%F0-标准化后的因变量数据
%T-自变量系统主成分得分n×rankE0矩阵
%h-用于建模或希望进行解释能力分析的主成分个数
%Rdx-各主成分对于某自变量的解释能力
%RdX-各主成分对于自变量组的解释能力
%RdXt-全部主成分对自变量组的累积解释能力
%Rdy-各主成分对于某因变量的解释能力
%RdY-各主成分对因变量组的解释能力
%RdYt-全部主成分对因变量组的累积解释能力

%--成分对自变量解释能力分析--
[nr,nx]=size(E0);
[nr,ny]=size(F0);
Rdx=zeros(nx,h);
t1=zeros(nr,1);
x1=zeros(nr,1);
for xj=1:nx
    for ti=1:h
        t1=T(:,ti);
        x1=E0(:,xj);
        cc=(corrcoef(t1,x1)).^2;
        Rdx(xj,ti)=cc(1,2);
    end
end
Rdx;
RdX=sum(Rdx)./nx;
RdXt=sum(RdX);
%--成分对因变量解释能力的分析--
Rdy=zeros(ny,h);
y1=zeros(nr,1);
for yj=1:ny
    for ti=1:h
        t1=T(:,ti);
        y1=F0(:,yj);
        rr=(corrcoef(t1,y1)).^2;
        Rdy(yj,ti)=rr(1,2);
    end
end
Rdy;
RdY=sum(Rdy)./ny;
RdYt=sum(RdY);


















