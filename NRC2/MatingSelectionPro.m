function MatingPool = MatingSelectionPro(Population,ArcPop,N)
% SelectedIndex = zeros(1,N);

  MatingPool = [];
if length(ArcPop)<N
    %SelectedIndex= TournamentSelection(2,N,-sum(max(0,Population.cons),2));
    SelectedIndex= TournamentSelection(2,N,DensityCal(Population.decs));
    MatingPool = Population(SelectedIndex);
else
    [FrontNo1,~] = NDSort(Population.objs,Population.cons,inf);
    Density1 =  DensityCal(Population.objs);
    SelectedIndex1 = TournamentSelection(2,N,FrontNo1,Density1);  
    [FrontNo2,~] = NDSort(-ArcPop.objs,inf);
    Density2 =  DensityCal(ArcPop.objs);
    SelectedIndex2 = TournamentSelection(2,N,FrontNo2,Density2);
   % SelectedIndex2 = TournamentSelection(2,N,sum(max(0,ArcPop.cons),2));
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
     [N,M] = size(PopObj);  
    Zmin       = min(PopObj,[],1);
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(max(PopObj),N,1)-repmat(Zmin,N,1)+1e-20);
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = Inf;
    DistanceSort = sort(Distance,2);
    Density = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));  
 end
