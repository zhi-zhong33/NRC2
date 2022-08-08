function MatingPool = MatingSelectionElite(Population,ArcPop,N)

 MatingPool = [];
 
if length(ArcPop)<N ||sum(sum(max(0,Population.cons),2)==0)<N
    PopAll = [Population,ArcPop];
    [~,b]=unique(PopAll.objs,'rows');
    PopAll = PopAll(1,b);
    SelectedIndex= TournamentSelection(2,N,DensityCal(PopAll.decs));
    MatingPool= PopAll(SelectedIndex);
else
    Density1 =  DensityCal(Population.objs);
    Density2 =  DensityCal(ArcPop.objs);
    index1 = randi(N,1,N);
    index2 = randi(N,1,N);
    for i = 1:N 
        PopObj1 = Population.objs;
        PopObj2 = ArcPop.objs; 
       
        CV = sum(max(0,ArcPop.cons),2)==0;
        
        if rand < 0.65
          if  PopObj1(index1(i),:)<PopObj1(index2(i),:)
            MatingPool = [MatingPool, Population(index1(i))];
          elseif  PopObj1(index1(i),:)>PopObj1(index2(i),:)
             MatingPool = [MatingPool, Population(index2(i))];
          elseif Density1(index1(i)) < Density1(index2(i))
             MatingPool = [MatingPool, Population(index1(i))];
          else
             MatingPool = [MatingPool, Population(index2(i))];
          end
        else  
            if min(PopObj2(index1(i),:)<PopObj2(index2(i),:))==1 || min(PopObj2(index1(i),:)>PopObj2(index2(i),:))==1
                if CV(index1(i),:)< CV(index2(i),:)
                 MatingPool = [MatingPool, ArcPop(index1(i))];
                else
                 MatingPool = [MatingPool, ArcPop(index2(i))];    
                end
            elseif Density2(index1(i)) < Density2(index2(i))
                MatingPool = [MatingPool, ArcPop(index1(i))];
            else
                 MatingPool = [MatingPool, ArcPop(index2(i))];   
            end
            
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
