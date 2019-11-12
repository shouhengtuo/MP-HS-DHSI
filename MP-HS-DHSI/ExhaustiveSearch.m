clear;

Dim = 100;
sample_num = 3000;
cacheSize = 10;
epi_dim = 2;
startIndex = 1; endIndex = 100; 


   

           folder = 'ExhaustiveMDR_resultData\';
           dataFile = strcat(folder,['ExhaustiveMultiPower', num2str(Dim),'_SNPs']);
%            dataFile = strcat(folder,['MultiPower5_', num2str(Dim),'_SNPs__errorRate_',num2str(pre_error_rate),'.xls']);
           dataFileTime = strcat(folder,['Runtime5_', num2str(Dim),'.xls']);
           dataFileFEs = strcat(folder,['FEs5_', num2str(Dim),'.xls']);
            A = {'TIME','RUN TIME','Power'};
        sheet = 1;
       xlRange = 'b1';
       xlswrite(dataFile,A,sheet,xlRange)

    for dataIndex =[49:60]%1:72 %
    %                 epi_dim = 2;
    %                 startIndex = 1; endIndex = 100;       
            switch(dataIndex)
%                      case 1
%                         epi_dim = 3;
%                         startIndex = 1; endIndex = 100;
%                         filepath = '.\threewayBests\'; filename = 'threewayBests'; 
%                     case 2
%                         epi_dim = 4;
%                         startIndex = 101; endIndex = 200; 
%                         filepath = '.\fourwayBests\';filename = 'fourwayBests';
%                     case 3
%                         epi_dim = 4;
%                         startIndex = 1101; endIndex = 1200;
%                         filepath = '.\fourwayNoLowBests\';    filename = 'fourwayNoLowBests';
%                     case 4
%                         epi_dim = 5;
%                         startIndex = 201; endIndex = 300;
%                         filepath = '.\fivewayBests\';  filename = 'fivewayBests';
                      %% 800 samples
            case 1  
                model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.05';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.05\400CASE_EDM-1\400CASE_EDM-1_';
               
            case 2
                model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.1';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.1\400CASE_EDM-1\400CASE_EDM-1_';
            case 3
                model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.2';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.2\400CASE_EDM-1\400CASE_EDM-1_';
            case 4
                model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.5';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.5\400CASE_EDM-1\400CASE_EDM-1_';
                
            case 5
               model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.05';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.05\400CASE_EDM-1\400CASE_EDM-1_';                
            case 6
                model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.1';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.1\400CASE_EDM-1\400CASE_EDM-1_';
            case 7
                model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.2';
                 filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.2\400CASE_EDM-1\400CASE_EDM-1_';
            case 8
                model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.5';
               filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.5\400CASE_EDM-1\400CASE_EDM-1_';
               
            case 9
                model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.05';
                 filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.05\400CASE_EDM-1\400CASE_EDM-1_';
            case 10
                model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.1';
                filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.1\400CASE_EDM-1\400CASE_EDM-1_';
            case 11
                model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.2';
                 filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.2\400CASE_EDM-1\400CASE_EDM-1_';
            case 12
                model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.5';
                 filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.5\400CASE_EDM-1\400CASE_EDM-1_';
           %% 4000 samples
            case 49
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.05';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.05\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 50
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.1';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.1\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 51
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.2';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.2\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 52
                    model='Model-1';parameter='H2=0.005,PD=0.1,MAF=0.5';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model1\h2=0.005,Pd=0.1,MAF=0.5\2000CASE_EDM-1\2000CASE_EDM-1_';

                case 53
                   model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.05';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.05\2000CASE_EDM-1\2000CASE_EDM-1_';                
                case 54
                    model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.1';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.1\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 55
                    model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.2';
                     filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.2\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 56
                    model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.5';
                   filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\Model2\h2=0.02,Pd=0.1,MAF=0.5\2000CASE_EDM-1\2000CASE_EDM-1_';

                case 57
                    model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.05';
                     filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.05\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 58
                    model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.1';
                    filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.1\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 59
                    model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.2';
                     filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.2\2000CASE_EDM-1\2000CASE_EDM-1_';
                case 60
                    model='Model-3';parameter='H2=0.02,PD=0.1,MAF=0.5';
                     filepath='O:\No.4 [NTFS 算法代码]\根目录\tsh msg\algorithm code0\MOFHSST\Mode data\model3\h2=0.02,PD=0.1,MAF=0.5\2000CASE_EDM-1\2000CASE_EDM-1_';   



                 
                        
                        
            end
         
          
            pvalue = 1/nchoosek(Dim,epi_dim);
            pvalue2 = 1e-4;
%%
               DATAS =[];% round(rand(sampleSize,Dim-100)*(3-0.00001) - 0.5 +0.000001);

            FdataNum = 0; 

            power = 0;
            powerSet = [];
            CX = [Dim-1,Dim];
            Evalutation_Times = [];
            TIME = [];
           

         for dataSetId = startIndex:endIndex
                 %% LOAD DATA for DME MODELS
                    if dataSetId<10
                        noId = strcat('00',num2str(dataSetId));
                     elseif dataSetId<100
                        noId = strcat('0',num2str(dataSetId));
                     else
                        noId = num2str(dataSetId);
                     end
                    
                     data = dlmread(strcat(filepath,noId,'.txt'),'\t',1,0);
                    
                % fprintf(' \n   dataId: %3d ：  ', dataSetId); 
%          %% load functional SNPs data  
%                 a = dlmread(strcat(filepath,strcat('best',num2str(dataSetId),'.txt')),'\t',1,0);
%         %% Generating nonfunctional SNP data according to HW
%                 AA = 0; Aa = 0; aa = 0;
%                 for i =1:1500%1:sample_num
%                     for j = 1:epi_dim
%                         if a(i,j) == 2
%                             aa = aa + 1;
%                         elseif a(i,j) == 1
%                             Aa = Aa + 1;
%                         else
%                             AA = AA + 1;
%                         end
%                     end
%                 end
%                 AA = 2*AA /(sample_num*epi_dim);
%                 Aa = 2*Aa /(sample_num*epi_dim);
%                 aa = 2*aa /(sample_num*epi_dim);
%                 
%                 b = zeros(sample_num,Dim-epi_dim);
%                 for i = 1 : sample_num
%                     for j = 1:Dim-epi_dim
%                         r = rand;
%                         if r <= aa
%                             b(i,j) = 2;
%                         elseif r <= aa+Aa
%                             b(i,j) = 1;
%                         else
%                             b(i,j) = 0;
%                         end
%                           % b(i,j) = fix(rand*3);
%                     end
%                 end
%                 data = [b,a];

 
                
               
     
                  tic       
%                             %% G-TEST :
%                                Gtest_Set = [];
%                                for i = 1:98
%                                    for j = i+1:99
%                                        for k = j+1:100
%                                               Gtest_Pvalue = Gtest_score(data(:,[i,j,k]),data(:,Dim+1));
%                                                if Gtest_Pvalue < pvalue
%                                                    Gtest_Set = [Gtest_Set; [i,j,k]];                          
%                                                end
%                                        end
%                                    end
%                                end           
                 %% G-TEST :
                               Gtest_Set = [];
                               for i = 1:99
                                   for j = i+1:100
                                       
                                              Gtest_Pvalue = Gtest_score(data(:,[i,j]),data(:,Dim+1));
                                               if Gtest_Pvalue < pvalue
                                                   Gtest_Set = [Gtest_Set; [i,j]]  ;                                                      
                                               end
                                               
                                                if [i,j] == CX 
                                                  
                                                    if Gtest_Pvalue < pvalue
                                                       powerSet = [powerSet,1];
                                                       power = power + 1;
                                                    else
                                                        powerSet = [powerSet,0];
                                                    end
                                                 end
                                       
                                   end
                               end    



         time = toc ;
         TIME = [TIME, time];       

         end
          fprintf("data %3d: detection power = %5.3f, mean time:%f\n",dataIndex, power, mean(TIME));
               FES = nchoosek(Dim,epi_dim);
        Results = [mean(TIME),FES,power/100];
          sheet = 1;
           xlRange = strcat('B', num2str(dataIndex+1), ': D', num2str(dataIndex+1)) ;
           xlswrite(dataFile,Results,sheet,xlRange)
    end
   
       

     