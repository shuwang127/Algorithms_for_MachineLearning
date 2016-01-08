%随机划分训练样本集和测试样本集
function [name_tr,name_test]=randQ(x)
j=1;k=1;
tr(250,252)=0;
test(250,252)=0;
for i=1:500
   if(rand>=0.5)
       if (j>250)
           for a=1:252
           test(k,a)=x(i,a);
           end
           k=k+1;
       else
          for a=1:252
          tr(j,a)=x(i,a);
          end
          j=j+1;
       end
   else
        if(k>250)
            for a=1:252
            tr(j,a)=x(i,a);
            end
            j=j+1;
        else
            for a=1:252
            test(k,a)=x(i,a);
            end
            k=k+1;
        end
    end 
end
name_tr=tr;
name_test=test;


       

