function [Feasible,Infeasible] = EnvironmentalSelectionNew(Population,N)
      
      % 去除重复解
      [~,b]=unique(Population.objs,'rows');
      if size(b,1)> N
           Population = Population(1,b);
      end
      % 把约束当作一个目标，分层
      [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);
      
      % 提取其中的可行解，并且获得可行解的编号
      FeasibleIndex = sum(max(0,Population.cons),2)==0;
      Feasible = Population(FeasibleIndex);
      FrontNoFeasible= FrontNo(FeasibleIndex);
      Feasible = FeasibleUpdate(Feasible,FrontNoFeasible,N);
      
      % 提取不可行解，并且获得相应编号
      InfeasibleIndex = sum(max(0,Population.cons),2)~=0;
      Infeasible = Population(InfeasibleIndex);
      FrontNoInfeasible= FrontNo(InfeasibleIndex);
      Infeasible = InfeasibleUpdate(Infeasible,FrontNoInfeasible,N);
end

function Population = FeasibleUpdate(Population,FrontNo,N)

if length(Population)>N
    % 给层排序，确定最大层
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N);
    
    % 确定前面的层
    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    
    % 确定最后一层并删除多余的个体
    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    if sum(Last)~= N-sum(Next)
    Population2 = Truncation(Population2,Population,N-sum(Next));
    end
    %确定最后的层
    Population = [Population1,Population2]; 
end
end

function Population = InfeasibleUpdate(Population,FrontNo,N)

if length(Population) > N
    %% 给层排序，确定最大层
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N+1);
    %% 确定前面的层
    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    
    %% 确定最后一层的个体和所需要的个体
    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    N2 = N- sum(Next);
    %% Boudary sorting 
    [FrontNo2,MaxFNo2] = NDSort(-Population2.objs,N2);
    Next2 = FrontNo2 < MaxFNo2;
     
    %% Advanced crowding distance sorting
    CrowdDis = AdvancedCrowdingDistance(Population2.objs,sum(max(0,Population2.cons),2),FrontNo2);
    % CrowdDis = DensityCal(Population2)
   
    %% Select the solutions in the last front based on their crowding distances
    Last2     = find(FrontNo2==MaxFNo2);
    [~,Rank] = sort(CrowdDis(Last2),'descend');
    Next2(Last2(Rank(1:N2-sum(Next2)))) = true;
    Population2 = Population2(Next2);
    
    %% Population for next generation
    Population = [Population1,Population2];
    
end

end


function Population = Truncation(Population,PopAll,N)
% Select part of the solutions by truncation

    %% Truncation
    Zmin       = min(PopAll.objs,[],1);
   % Zmax = max(Population.objs,[],1);
    PopObjTemp = Population.objs;
    PopObj = (PopObjTemp -repmat(Zmin,length(Population),1))./(repmat(max(PopAll.objs),length(Population),1)-repmat(Zmin,length(Population),1));
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
    
      Population = Population(Del==0);
end


function CrowdDis = AdvancedCrowdingDistance(PopObj,CV,FrontNo)
% Calculate the crowding distance of each solution front by front
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);
        Fmin  = min(PopObj(Front,:),[],1);
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
            CrowdDis(Front(Rank(1)))   = inf;
            CrowdDis(Front(Rank(end))) = inf;
            for j = 2 : length(Front)-1
              CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
             % CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
 %  CrowdDis = CrowdDis./(1+CV);
end

