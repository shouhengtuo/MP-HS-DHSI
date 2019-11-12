function [P_value] = Gtest_score(snp_com,state)
%compute the chi-score
P_value = 1;
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end

[xrow,xcol] = size(snp_com);


[Data,idx,cid]=unique(snp_com,'rows');
[lrow,~]=size(Data);


%% 

F=zeros(3,lrow+1);  %% 第一行存放 case 数， 第二行存放 control 数， 第三行存放 总数
E=zeros(2,lrow);  %% 期望数

for i=1:xrow   %% 统计每个基因型组合出现的次数
   
   
       if state(i)==0 %% case        
            F(1,cid(i))=F(1,cid(i))+1;        
       else 
            F(2,cid(i))=F(2,cid(i))+1;  
       end  
    
end



F(3,1:lrow)=sum(F(1:2,1:lrow),1);
F(:,lrow+1)=sum(F(:,1:lrow),2);


G=0;
Degree=(2-1)*(lrow-1);

for i=1:2
    for j=1:lrow
        O=F(i,j)/xrow;
        E=( ( F(i,lrow+1) * F(3,j) )/xrow )/xrow;
       
      if F(3,j)>1
          if O>0
              G=G+(xrow*O)*log(O/E);
              %G=G+O*log(O/E);
          end
      
      elseif Degree>1
         Degree=Degree-0.5; 
      end
    end    
end

G=2*G;

P_value=1-chi2cdf(G,Degree);


