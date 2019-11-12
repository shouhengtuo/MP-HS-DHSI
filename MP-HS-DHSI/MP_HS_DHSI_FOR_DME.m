clear;

Dim = 100;
sample_num = 4000;
cacheSize = 20;
epi_dim = 2;
startIndex = 1; endIndex = 100; 

pre_error_rate = 0.45;

    %% harmony search algorithm parameters setting
          pvalue = 0.05/nchoosek(Dim,epi_dim);

          
           HMS = 100;
           CandidateSize =5;

           
    CX =[99 100];

           folder = 'DME_result\';
    %        dataFile = strcat(folder,['MultiPower', num2str(Dim),'_SNPs']);
           dataFile = strcat(folder,['MultiPower', num2str(Dim),'_SNPs__errorRate_',num2str(pre_error_rate),'.xls']);
            
           dataFileTime = strcat(folder,['Runtime_', num2str(Dim),'.xls']);
           dataFileFEs = strcat(folder,['FEs_', num2str(Dim),'.xls']);
            A = {'1st Power','2nd Power','3rd Power', 'K2_Power', 'Gtest_Power','LR_Power','Gini_Power', 'TPR','SPC', 'PPV', 'ACC','FDR','F1', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime','Power4', 'TPR2','SPC2', 'PPV2', 'ACC2','FDR2','F1-2'};
        sheet = 1;
       xlRange = 'b1';
       xlswrite(dataFile,A,sheet,xlRange)


    for dataIndex =1:4
    %                 epi_dim = 2;
    %                 startIndex = 1; endIndex = 100;       
            switch(dataIndex)

           
    %% 4000 samples
                 case 1
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.05';
                    filepath='.\DME_DATA\DMEmodel1_11_EDM-1\DMEmodel1_11_EDM-1_';
                case 2
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.1';
                    filepath='.\DME_DATA\DMEmodel1_12_EDM-1\DMEmodel1_12_EDM-1_';
                 case 3
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.05';
                    filepath='.\DME_DATA\DMEmodel1_13_EDM-1\DMEmodel1_13_EDM-1_';
                case 4
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.1';
                    filepath='.\DME_DATA\DMEmodel1_14_EDM-1\DMEmodel1_14_EDM-1_';
            end
             max_iter =3000; % min([Dim * 100,100000]);
           maxIterForLocalSearch = 1000;
            
         

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
                 if dataSetId<10
                    noId = strcat('00',num2str(dataSetId));
                 elseif dataSetId<100
                    noId = strcat('0',num2str(dataSetId));
                 else
                    noId = num2str(dataSetId);
                 end
                 %% load data
                 % data format
                 % line: snp1, snp2, ..., snpN, label
                data = dlmread(strcat(filepath,noId,'.txt'),'\t',1,0);
               

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
                     fprintf('\nIndex:%d, dataSet %3d:search time(%f)|FEs=%d,   %d/%d    successfully(^V^)    ',dataIndex,dataSetId, runtime,Nc,power1,dataSetId-startIndex+1);
                else 
                    fprintf('\nIndex:%d, dataSet %3d: search time(%f)|FEs=%d   %d/%d   search failed!!!!!! ! ',dataIndex,dataSetId,runtime,Nc,power1,dataSetId-startIndex+1);
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
      
       G_set = [];
                for i = 1:G_size
                            [CERR_RATE,PERR_RATE,CA,BA] = MDR2(data(:,Gtest_Set(i,:)),data(:,Dim+1));
         
                           if  PERR_RATE <pre_error_rate ||  BA > 1- pre_error_rate 
                                G_set = [G_set;Gtest_Set(i,1:epi_dim)]; 
                           end
%                            if Gtest_Set(i,1:epi_dim) == CX
%                                fprintf('**[%5.3f,%5.3f,%5.3f,%5.3f]   ',  CERR_RATE,PERR_RATE,CA,BA);                               
%                            end
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