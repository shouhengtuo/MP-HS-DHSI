
function minDist=getMinDistance2(cc,Centers,dim_epi)
       [m,n] = size(Centers);
      X = reshape(Centers,m*n,1);
      CF = tabulate(X); % 统计重复元素出现的次数
      [cx,idx] = sort(CF(:,2),'descend');
      aX = CF(idx(1:dim_epi),1);
     
      center = intersect(cc,aX');
      minDist = length(center);    
 end