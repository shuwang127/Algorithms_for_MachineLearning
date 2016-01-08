function cr=plsutcor(U,T) 
%PLSUTCOR绘制t1/u1图，并给出二者的相关系数
%cr=plsutcor(U,T)
%U-因变量提取的成分
%T-自变量提取的成分
%cr-自变量与因变量的相关系数

u1=U(:,1);
t1=T(:,1);
ut=[u1,t1];
cr=corrcoef(ut)
plot(t1,u1,'o')
lsline
title('t1/u1图')
xlabel('t1')
ylabel('u1')
