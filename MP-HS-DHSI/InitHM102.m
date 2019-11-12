function [X,HM,Fit]=InitHM102(data,HMS,dim_epi,Center)
    X=zeros(HMS,dim_epi);
    SNPs=length(data(1,:))-1;
    snp=[];
    for i=1:HMS

        snp(1)=ceil(rand*SNPs);
        for j=2:dim_epi
          snp(j)=ceil(rand*SNPs);  

          while ismember(snp(j),snp(1:j-1)) 
             snp(j)=ceil(rand*SNPs);        
          end
        end
        Temp=snp;
        snp=sort(snp);
        
        c= 0;
        while ismember(snp,X,'rows') || getMinDistance2(snp,Center(:,1:dim_epi),dim_epi)<=dim_epi-1
            j=ceil(rand*dim_epi);
            ss=ceil(rand*SNPs);  
            while( ismember(ss,snp))
                ss=ceil(rand*SNPs);  
            end
            snp(j) = ss;
            Temp=snp;
            snp=sort(snp);
            c = c + 1;
            if  c > dim_epi*5
                break;
            end
        end
        
        
        
        HM(i,:)=Temp;
        X(i,:)=snp;
%         size(data)
%         X(i,:)
%         SNPs
        [Fit(i,1),Fit(i,2),Fit(i,3),Fit(i,4)] = multi_criteriaEvaluationFuns3(data(:,X(i,:)),data(:,end)); 
        
        snp=[];
        
%            if  c > dim_epi*10
%                 break;
%             end
    end
    
end