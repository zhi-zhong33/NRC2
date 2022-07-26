function ArcPop = UpdateArc(Population,N)

    Population = Population(sum(max(0,Population.cons),2)==0);
    
     
     [~,b]=unique(Population.objs,'rows');
     Population = Population(1,b);
           
             
    [FrontNo,~] = NDSort(Population.objs,1);
    Temp = FrontNo==1;
    Population = Population(Temp==1);
    if length(Population)<=N 
      ArcPop = Population;
    else 
       ArcPop = Truncation(Population,N);    
    end 
end

% end

function ArcPop = Truncation(Population,N)
% Select part of the solutions by truncation

    %% Truncation
    Zmin       = min(Population.objs,[],1);
   % Zmax = max(Population.objs,[],1);
    PopObj = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(Population.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1)+1e-10)+1e-10;
   % PopObj = (PopObj(Last,:)); 
    
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < size(PopObj,1)-N
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
    
      ArcPop = Population(Del==0);
end


