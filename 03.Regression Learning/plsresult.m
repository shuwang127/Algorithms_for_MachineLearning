function result = plsresult( W,RdY,RdYt,h )

[nx,wk]=size(W);
result=zeros(1,nx);
for j=1:nx
    for hh=1:h
        Whj=W(j,hh);
        tvip(hh)=RdY(hh)*Whj.^2;
    end
    S_tvip=sum(tvip);
    result(j)=sqrt((nx/RdYt)*S_tvip);
end
bar(result,'c');
title('变量投影重要性result图')