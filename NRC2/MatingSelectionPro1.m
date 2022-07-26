function MatingPool = MatingSelectionPro1(Population,ArcPop,N)
% SelectedIndex = zeros(1,N);

  MatingPool = [];
   
%if length(ArcPop)<N || length(Population)<N
if length(ArcPop)<N ||sum(sum(max(0,Population.cons),2)==0)<N
    PopAll = [Population,ArcPop];
    [~,b]=unique(PopAll.objs,'rows');
    PopAll = PopAll(1,b);
    SelectedIndex= TournamentSelection(2,N,DensityCal(PopAll.decs));
    MatingPool= PopAll(SelectedIndex);
else
    Density1 =  DensityCal(Population.objs);
    SelectedIndex1 = TournamentSelection(2,N,Density1);  
    Density2 =  DensityCal(ArcPop.objs);
    SelectedIndex2 = TournamentSelection(2,N,Density2);
    for i = 1:N
        if rand < 0.8
             MatingPool = [MatingPool,Population(SelectedIndex1(i))];
        else
             MatingPool = [MatingPool,ArcPop(SelectedIndex2(i))];
        end
    end
    
end
end  
 


 function Density = DensityCal(PopObj)
     [N,~] = size(PopObj);  
    Zmin       = min(PopObj,[],1);
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(max(PopObj),N,1)-repmat(Zmin,N,1)+1e-20);
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = Inf;
    DistanceSort = sort(Distance,2);
    Density = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));  
 end
