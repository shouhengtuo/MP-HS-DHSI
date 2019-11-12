function [CERR_RATE,PERR_RATE,CA,BA] = MDR2(SNP_COM,state)
 % SNP_COM 基因型 ： 0 1 2 
% state   表型 1 case ； 0 control

%%
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end

[xrow,xcol] = size(SNP_COM);
[Data,idx,cid]=unique(SNP_COM,'rows');
[lrow,~]=size(Data);
sample=zeros(lrow,1);
disease=sample;
control=sample;

for i=1:xrow   %% 统计每个基因型组合出现的次数
   if state(i) == 1
       disease(cid(i)) = disease(cid(i)) + 1;
   else
       control(cid(i)) = control(cid(i)) + 1;
   end
end
sample = disease + control;
ProCaseToControl = (disease+1)./(control + 1);
% threshold = sum(ProCaseToControl.* (sample/xrow));
 sample_num = xrow;
caseNum = sum(state);
controlNum = sample_num - caseNum;
threshold = caseNum/controlNum;
%%

 epi_dim  = xcol;


caseId = find(state==1);
controlId = find(state == 0);
CaseSample = SNP_COM(caseId,:);
ControlSample = SNP_COM(controlId,:);


for c = 1:10 % 10倍交叉验证
    caseTrainNum = fix(caseNum*0.9);
    controlTrainNum = fix(controlNum*0.9);
    a = randperm(caseNum);
    caseTrain = CaseSample(a(1:caseTrainNum),:);
    caseTest = CaseSample(a(caseTrainNum+1:caseNum),:);
    a = randperm(controlNum);
    controlTrain = ControlSample(a(1:controlTrainNum),:);
    controlTest = ControlSample(a(controlTrainNum+1:controlNum),:);
    

    
        caseTrainAcc = zeros(3,3);
        controlTrainAcc = caseTrainAcc;
        
    for i = 1:caseTrainNum
        A = caseTrain(i,:) + 1;
        caseTrainAcc(A(1),A(2)) = caseTrainAcc(A(1),A(2)) + 1;
        A = controlTrain(i,:) + 1;
        controlTrainAcc(A(1),A(2)) = controlTrainAcc(A(1),A(2)) + 1;
    end
 caseTestAcc = zeros(3,3);
controlTestAcc = caseTestAcc;
    for i = 1:caseNum - caseTrainNum
        A = caseTest(i,:) + 1;
        caseTestAcc(A(1),A(2)) = caseTestAcc(A(1),A(2)) + 1;
        A = controlTest(i,:) + 1;
        controlTestAcc(A(1),A(2)) = controlTestAcc(A(1),A(2)) + 1;
    end
    
    riskTrainM = caseTrainAcc./controlTrainAcc;
%     %%
%     threshold = mean(riskTrainM);
    riskTestM = caseTestAcc./controlTestAcc;
    [m,n] = size(riskTrainM);
    %% 分类错误率估计
    C_TP = 0;  C_TN = 0; C_FN = 0; C_FP = 0;
    %% 预测错误率估计
    P_TP = 0;  P_TN = 0; P_FN = 0; P_FP = 0;
    
    HRcase = 0; HRcontrol = 0;
    LRcase = 0; LRcontrol = 0;
    
    for i = 1:m
        for j = 1:n
            if riskTrainM(i,j) >= threshold
                HRcase = HRcase + caseTrainAcc(i,j);
                HRcontrol = HRcontrol + controlTrainAcc(i,j);
                C_TP = C_TP + caseTrainAcc(i,j);C_FP = C_FP + controlTrainAcc(i,j);
                 if riskTestM(i,j) >= threshold
                     P_TP = P_TP + caseTestAcc(i,j);  P_FP = P_FP + controlTestAcc(i,j);
                 else
                     P_FN = P_FN + caseTestAcc(i,j); P_TN = P_TN + controlTestAcc(i,j);
                 end
            else
                LRcase = LRcase + caseTrainAcc(i,j);
                LRcontrol = LRcontrol + controlTrainAcc(i,j);
                C_TN = C_TN + controlTrainAcc(i,j);C_FN = C_FN + caseTrainAcc(i,j);
                if riskTestM(i,j) >= threshold
                     P_FP = P_FP + controlTestAcc(i,j);P_TP = P_TP + caseTestAcc(i,j); 
                else
                   P_TN = P_TN + controlTestAcc(i,j); P_FN = P_FN + caseTestAcc(i,j);
                end
            end
        end
    end
    ConTable = [HRcase HRcontrol; LRcase LRcontrol];
    ConTable(3,:) = sum(ConTable);
    ConTable(:,3) = sum(ConTable,2);
    ExpTable(1,1:2) = ConTable(3,1:2)*ConTable(1,3)/ConTable(3,3);
    ExpTable(2,1:2) = ConTable(3,1:2)*ConTable(2,3)/ConTable(3,3);
    
    LR = 2*sum(sum(ConTable(1:2,1:2).*log10(ConTable(1:2,1:2)./ExpTable)));
    X2 = sum(sum((((ConTable(1:2,1:2)./ExpTable).^2)./ExpTable)));
%     1-chi2cdf(2*X2,1)
     CC_TP(c) = C_TP/caseTrainNum;  
    CC_TN(c) = C_TN/controlTrainNum;
    CC_FN(c) = C_FN/caseTrainNum;
    CC_FP(c) = C_FP /controlTrainNum;
    CA(c) = (CC_TP + CC_TN)/(CC_TP + CC_TN + CC_FN + CC_FP);   %% 分类精度
    
    PP_TP(c) = P_TP/(caseNum - caseTrainNum);  
    PP_TN(c) = P_TN/(controlNum - controlTrainNum);
    PP_FN(c) = P_FN/(caseNum - caseTrainNum);
    PP_FP(c) = P_FP /(controlNum - controlTrainNum);
    BA(c) = 0.5*(PP_TP/(PP_TP+PP_FN) + PP_TN/(PP_TN+PP_FP)); %% 预测精度
    

end
mean([CC_TP; CC_TN; CC_FN; CC_FP],2);
CERR_RATE = 0.5*(CC_FN/(CC_TP+CC_FN) + CC_FP/(CC_FP+CC_TN));  %% classification error rate
CA = mean(CA); %% 分类精度
mean([PP_TP; PP_TN; PP_FN; PP_FP],2);
PERR_RATE = 0.5*(PP_FN/(PP_TP+PP_FN) + PP_FP/(PP_FP+PP_TN));  %% prediction error rate 
BA= mean(BA); %% 


   % score = mean(errorRate)*100;
%     rid = randperm(sample_num);
%     trainX = SNP_COM(rid(1:train_num),:);
%     trainY = state(rid(1:train_num));
%     testX = SNP_COM(rid(train_num + 1:sample_num),:);
%     testY = state(rid(train_num + 1:sample_num));
    
%% training




end