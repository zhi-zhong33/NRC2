function DensityValue = Density(Population)
    % DensityValue = [];
     Zmin       = min(Population.objs,[],1);
     Zmax = max(Population.objs,[],1);
    PopObj = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(Zmax,length(Population.objs),1)-repmat(Zmin,length(Population.objs),1));
   % PopObj = (PopObj(Last,:)); 
   % PopObj = real(PopObj);
   % PopObj
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
   % DistanceSort = sort(Distance,2);
   % DensityValue = DistanceSort(:,floor(sqrt(length(Distance)))+1); 
    DensityValue    = min(Distance,[],2);
end