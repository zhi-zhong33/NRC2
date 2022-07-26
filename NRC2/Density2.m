function DensityValue = Density2(Population)
 
 
       Zmin       = min(Population.objs,[],1);
 
    PopObj= (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(Population.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1));
    Cosine   = 1 - pdist2(PopObj,PopObj,'cosine');
    Cosine   = Cosine.*(1-eye(size(PopObj,1)));
    %DensityValue    = sort(Cosine2,2,'descend');
    DensityValue    = -max(Cosine,[],2);
end