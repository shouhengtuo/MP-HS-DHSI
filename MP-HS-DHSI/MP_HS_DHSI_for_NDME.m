clear;

Dim = 100;
sample_num = 3000;
cacheSize = 10;
epi_dim = 2;
startIndex = 1; endIndex = 100; 
pre_error_rate = 0.45;


    %% harmony search algorithm parameters setting          
          
           HMS = 100 % population size 
           CandidateSize =5;

           permutation_times = 100;
   

           folder = 'resultData\';
    %        dataFile = strcat(folder,['MultiPower', num2str(Dim),'_SNPs']);
           dataFile = strcat(folder,['MultiPower_', num2str(Dim),'_SNPs__errorRate_',num2str(pre_error_rate),'.xls']);
           dataFileTime = strcat(folder,['Runtime_', num2str(Dim),'.xls']);
           dataFileFEs = strcat(folder,['FEs_', num2str(Dim),'.xls']);
            A = {'1st Power','2nd Power','3rd Power', 'K2_Power', 'JS_Power','LR_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC','FDR','F1', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime','Power4', 'TPR2','SPC2', 'PPV2', 'ACC2','FDR2','F1-2'};
        sheet = 1;
       xlRange = 'b1';
       xlswrite(dataFile,A,sheet,xlRange)


    for dataIndex =[1:4]%1:72 %
    %                 epi_dim = 2;
    %                 startIndex = 1; endIndex = 100;       
            switch(dataIndex)
                     case 1
                        epi_dim = 3;
                        startIndex = 1; endIndex = 100;
                        filepath = '.\threewayBests\'; filename = 'threewayBests'; 
                    case 2
                        epi_dim = 4;
                        startIndex = 101; endIndex = 200; 
                        filepath = '.\fourwayBests\';filename = 'fourwayBests';
                    case 3
                        epi_dim = 4;
                        startIndex = 1101; endIndex = 1200;
                        filepath = '.\fourwayNoLowBests\';    filename = 'fourwayNoLowBests';
                    case 4
                        epi_dim = 5;
                        startIndex = 201; endIndex = 300;
                        filepath = '.\fivewayBests\';  filename = 'fivewayBests';
                    
                 
                        
                        
            end
         
            CX = Dim-epi_dim+1:Dim;
            max_iter = 3000 * epi_dim^3;

             maxIterForLocalSearch = max(2000,max_iter/3);% number of iteration for local search
             pvalue = 0.01/nchoosek(Dim,epi_dim);
             pvalue2 = 1e-4;
%%
               DATAS =[];% round(rand(sampleSize,Dim-100)*(3-0.00001) - 0.5 +0.000001);

            FdataNum = 0; 

            Ac = 0;
            TP = 0;
            FP = 0;
            TN = 0;
            FN = 0;

            TP2 = 0;
            FP2 = 0;
            TN2 = 0;
            FN2 = 0;
             power1 = 0;  
             power2 = 0; 
             power3 = 0;
             power4 = 0;
             K2_power = 0;
             JS_Power = 0;
             LR_power = 0;
             JE_Power = 0;
             JE_power = 0;
            Evalutation_Times = [];
            TIME = [];
             succ = 0; 
            succEvalutation_Times = [];
            succTIME = [];

         for dataSetId = startIndex:endIndex
                 %% LOAD DATA
    
                % fprintf(' \n   dataId: %3d ：  ', dataSetId); 
                %% load functional data
                a = dlmread(strcat(filepath,strcat('best',num2str(dataSetId),'.txt')),'\t',1,0);
                %% generate non-functiondata based HWD
                AA = 0; Aa = 0; aa = 0;
                for i =1:1500%1:sample_num
                    for j = 1:epi_dim
                        if a(i,j) == 2
                            aa = aa + 1;
                        elseif a(i,j) == 1
                            Aa = Aa + 1;
                        else
                            AA = AA + 1;
                        end
                    end
                end
                AA = 2*AA /(sample_num*epi_dim);
                Aa = 2*Aa /(sample_num*epi_dim);
                aa = 2*aa /(sample_num*epi_dim);
                
                b = zeros(sample_num,Dim-epi_dim);
                for i = 1 : sample_num
                    for j = 1:Dim-epi_dim
                        r = rand;
                        if r <= aa
                            b(i,j) = 2;
                        elseif r <= aa+Aa
                            b(i,j) = 1;
                        else
                            b(i,j) = 0;
                        end
                          % b(i,j) = fix(rand*3);
                    end
                end
                data = [b,a];


                
                 %% Search k-snp loci using Harmony search algorithm        
                 [Candidate,canSize,Nc,runtime,flag] = HS_2019_multiCRITERIA5(data,epi_dim,HMS,max_iter,maxIterForLocalSearch,CandidateSize,CX);
                 Evalutation_Times = [Evalutation_Times , Nc];
                 TIME = [TIME  runtime];
                 if fix(flag/1000) == 1
                     K2_power = K2_power + 1;
                 end

                 if fix(mod(flag,1000)/100) == 1                
                     JS_Power = JS_Power + 1;              
                 end
                 if fix(mod(flag,100)/10) == 1
                     LR_power = LR_power + 1;             
                 end
                 if fix(mod(flag,10)) == 1
                      JE_Power = JE_Power + 1;         
                 end
             % end harmony search  ****************************************
                  [SNP_COM1, bestId1] = sort(Candidate(:,epi_dim+1));
                  [SNP_COM2, bestId2] = sort(Candidate(:,epi_dim+2));
                  [SNP_COM3, bestId3] = sort(Candidate(:,epi_dim+3));


               %% 1st stage power  
                if flag > 0
                     power1 = power1 + 1;
                     succ = succ + 1;
                     succEvalutation_Times = [succEvalutation_Times, Nc];
                     succTIME = [succTIME  runtime];
                     fprintf('\nIndex:%d, dataSet %3d:search time(%f)|FEs=%d,   %d/%d    success(^V^)    ',dataIndex,dataSetId, runtime,Nc,power1,dataSetId-startIndex+1);
                else 
                    fprintf('\nIndex:%d, dataSet %3d: search time(%f)|FEs=%d   %d/%d   !!!!!! ! ',dataIndex,dataSetId,runtime,Nc,power1,dataSetId-startIndex+1);
                end

    %% G-TEST :
       Gtest_Set = [];
       for i = 1:canSize
           Gtest_Pvalue = Gtest_score(data(:,Candidate(i,1:epi_dim)),data(:,Dim+1));
           if Gtest_Pvalue < pvalue
               Gtest_Set = [Gtest_Set; Candidate(i,1:epi_dim)];
              if Candidate(i,1:epi_dim) == CX
                   power2 = power2 + 1;
              end
           end
       end           

    %% 3nd stage : MDR 
       if isempty(Gtest_Set)
           G_size = 0;
       else
          G_size = length(Gtest_Set(:,1));
       end
        fprintf('Gtest-set size = %d ,',  G_size);   
       G_set = [];
                for i = 1:G_size
                            [CERR_RATE,PERR_RATE,CA,BA] = MDR2(data(:,Gtest_Set(i,:)),data(:,Dim+1));
         
                           if  PERR_RATE <pre_error_rate ||  BA > 1- pre_error_rate 
                                G_set = [G_set;Gtest_Set(i,1:epi_dim)]; 
                           end
                           if Gtest_Set(i,1:epi_dim) == CX
                               fprintf('**[%5.3f,%5.3f,%5.3f,%5.3f]   ',  CERR_RATE,PERR_RATE,CA,BA);                               
                           end
                end

               if ~isempty(G_set) 
                   if ismember(CX,G_set,'rows')
                       power3 = power3 + 1;
                       TP = TP + 1;
                       TN = TN + 1;
                   else
                       FN = FN + 1;
                       FP = FP + 1;
                   end
                  G_setSize = length(G_set(:,1));

               else
                    G_setSize = 0;  
                    FP = FP + 1;
                    TN = TN + 1;
               end
                fprintf('MDR-set size = %d \n ',  G_setSize);

         end

        Power1 = power1/100;
       Power2 = power2/100;
       Power3 = power3/100;
       K2_Power = K2_power/100;
       JS_Power = JS_Power/100;
       LR_Power = LR_power/100;
       JE_Power = JE_Power/100;
%        JE_Power = JE_power/100;
       E_Times = median( Evalutation_Times);
       RunTime = median( TIME);
       succE_Times = median (succEvalutation_Times);
       succRunTime = median( succTIME);
    %    真阳性率(True Positive Rate，TPR)，灵敏度(Sensitivity)，召回率(Recall)：
    % Sensitivity=Recall=TPR

      TPR = TP/(TP + FN + 1e-10);
      SPC = TN/(FP + TN + 1e-10);
      PPV = TP/(TP + FP+ 1e-10);
      ACC = (TP + TN)/(TP + TN + FN + FP);
      FDR = 1 - PPV;
      F1 = 2*TP/(2*TP + FP + FN);
      
        Results=[Power1,Power2,Power3, K2_Power,JS_Power, LR_Power,JE_Power, TPR,SPC, PPV, ACC, FDR,F1, E_Times, RunTime, succE_Times, succRunTime];


       % A = {'Power1','Power2','Power3', 'K2_Power', 'LR_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime'};
        sheet = 1;
       xlRange = strcat('B', num2str(dataIndex+1)) ;
       xlswrite(dataFile,Results,sheet,xlRange)
       

       xlRange = strcat('B', num2str(dataIndex+1)) ;
       xlswrite(dataFileTime,TIME,sheet,xlRange);
       xlswrite(dataFileFEs,Evalutation_Times,sheet,xlRange);
    end
