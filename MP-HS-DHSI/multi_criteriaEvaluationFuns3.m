% function [K2Score,Gtest,LR,Gini_Score,JE_Score,ME] = multi_criteriaEvaluationFuns(snp_com,state)
 function [K2Score,JS,LR,JE_Score] = multi_criteriaEvaluationFuns3(snp_com,state)
% function [K2Score,chiT,LR,JS] = multi_criteriaEvaluationFuns(snp_com,state)
%%
% 
%%
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end

[xrow,xcol] = size(snp_com);
[Data,idx,cid]=unique(snp_com,'rows');
[lrow,~]=size(Data);

Hcase = zeros(3,xcol);
Hcontrol = zeros(3,xcol);
sample=zeros(lrow,1);
disease=sample;
control=sample;

for i=1:xrow   %% 统计每个基因型组合出现的次数
   if state(i) == 1
       disease(cid(i)) = disease(cid(i)) + 1;
         for j = 1:xcol
              if snp_com(i,j) == 0
                  Hcase(1,j) = Hcase(1,j) + 1;
              elseif snp_com(i,j) == 1
                  Hcase(2,j) = Hcase(2,j) + 1; 
              else
                  Hcase(3,j) = Hcase(3,j) + 1;
              end
         end
   else
       control(cid(i)) = control(cid(i)) + 1;
        for j = 1:xcol
              if snp_com(i,j) == 0
                  Hcontrol(1,j) = Hcontrol(1,j) + 1;
              elseif snp_com(i,j) == 1
                  Hcontrol(2,j) = Hcontrol(2,j) + 1; 
              else
                  Hcontrol(3,j) = Hcontrol(3,j) + 1;
              end
         end
   end
  
end

A = sum((Hcase - Hcontrol).^2);
sample = disease + control;
%[disease,control]
SDC = sqrt(sum(A))/sum(abs(disease-control));
 %% G-test
     F = [disease';control'];
     F(3,1:lrow)=sum(F(1:2,1:lrow),1);
     F(:,lrow+1)=sum(F(:,1:lrow),2);
     
     G=0;
     chi = 0;
    % LR = 0;
     Degree=(2-1)*(lrow-1);
    for i=1:2
        for j=1:lrow
            O = F(i,j)/xrow;
            E = ( ( F(i,lrow+1) * F(3,j) )/xrow )/xrow;

              if F(3,j)>1
                   if O>0
                      G = G+(xrow * O) * log(O/E);
                     % LR = LR + O * log(O/E);
                      chi = chi + (abs(O - E))^2 / E;
                   end

              elseif Degree>1
                 Degree = Degree-0.5; 
              end
        end    
    end
chi = 2 * chi;
% chiT = 1 / chi;
G=2*G;
%Gtest = 1/G;
LR = 1 / G;
Gtest = 1 - chi2cdf(G,Degree);
%1/LR;
%GtestP_value=1-chi2cdf(G,Degree);

% [disease';control';sample']
%% K2 score 
    K2score=0;

    for i=1:lrow
        if sample(i)>0
            y=My_factorial(sample(i)+1);
            r=My_factorial(disease(i))+My_factorial(control(i));
            K2score=K2score+(r-y);
        end
    end
   K2Score =abs(K2score);



%% GINI score
 

%         sCase = sum(disease);
%         Pcase = disease./sCase;%
%         Pcontrol = control./(xrow - sCase);%
%         P = sample / xrow;
% 
%          Gini_Score = sCase/xrow * (1 - sum(Pcase.^2)) / ((1-sCase/xrow) * (1-sum(Pcontrol.^2)));

        
  
 %% 互熵 mutual entropy
%      Psample = sample / xrow;
%      PCC = [disease; control]/xrow;
%     
%      PCC = PCC(PCC>0);
%      s1 = sCase/xrow; s2 = 1 - s1;
%      MeY = - (s1 .* log2(s1) + s2 .* log2(s2));
%    
%      MeX = - sum(Psample.*log2(Psample));
%      MeXY = - sum ( PCC .* log2(PCC));
%      ME = 1 / ( MeY + MeX - MeXY);


%% Joint Entropy of disease genotype combinantion 
        Psample = sample / xrow;
        Pdisease = disease / sum(disease);
        Pcontrol = control / sum(control);
        JE = 0;
        JC = 0;
        JCC = 0;
        for i = 1:lrow
           if Pdisease(i)>0
              JE = JE + Pdisease(i).*log2(Pdisease(i));
           end
            if Pcontrol(i)>0
              JC = JC + Pcontrol(i).*log2(Pcontrol(i));
            end
%             JCC = JCC + Pdisease(i).*log2(Pdisease(i)) / (Pcontrol(i).*log2(Pcontrol(i)));
            
        end
%          JE_Score = JointEntropyScore * JC;
%          JE_Score = -JointEntropyScore/JE_Score;
      
        %  JE_Score = -JE/ JCC^2  ;
%          JE_Score =  -JE / JC^2;
         % JE_Score = - 1/JC;
%           JE_Score = -JE / ((JE - JC)^2); % - 1/JC;

            
%            JE_Score = SDC / (JE + JC)^2 ;% (JE + JC)^2  is quadratic sum of the joint entropy of SNP combintion
          JE_Score = SDC / JC^2 ;
     %     JE_Score = SDC;
   
 
          
          
          %% Jensen-Shannon
 JS1 = 0;
 JS2 = 0;
 JF = 0; %Jeffrey
     for i = 1:lrow
           if Pdisease(i)>0
              JS1 = JS1 + Pdisease(i).*log2(2*Pdisease(i) / (Pdisease(i) + Pcontrol(i)) );
           end
            if Pcontrol(i)>0
              JS2 = JS2 + Pcontrol(i).*log2(2*Pcontrol(i) / (Pdisease(i) + Pcontrol(i)) );
            end
           
%             if Pdisease(i)+Pcontrol(i) > 0
%                 JF = JF + ((Pdisease(i)-Pcontrol(i))^2) / Pcontrol(i);
%             end
     end
      %JS =1- 0.5*(JS1 + JS2) ;
      JS = 1 /(JS1 + JS2);
      
%% Wasserstein Distance
%    WS = 1 / sum((Pdisease-Pcontrol).^2);
 

    


%% f is function used to calculate log form factorial
    function f=My_factorial(e)
        f=0;
        if e>0
            for o=1:e
                f=f+log(o);
            end
        end
    end
end


