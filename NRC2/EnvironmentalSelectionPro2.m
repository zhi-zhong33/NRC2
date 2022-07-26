function [Feasible,Infeasible] = EnvironmentalSelectionPro2(Population,N)
      
%       % 去除重复解
%       [~,b]=unique(Population.objs,'rows');
%       Population = Population(1,b);
      
      [FrontNo,MaxNo] = NDSort([Population.objs,sum(max(0,Population.cons),2)],2*N);

      % 提取其中的可行解，并且获得可行解的编号
      FeasibleIndex = sum(max(0,Population.cons),2)==0;
      FeasibleFrontNo = FrontNo(FeasibleIndex);
      NextFeasibleIndex = FeasibleFrontNo<=MaxNo;
      FeasibleIndex1 = FeasibleIndex(NextFeasibleIndex);
      Feasible = Population(FeasibleIndex1);
      FeasibleFrontNo = FrontNo(FeasibleIndex1); 
      Feasible = FeasibleUpdate(Feasible,FeasibleFrontNo,N);
      
          
    
      % 提取不可行解，并且获得相应编号
      InfeasibleIndex = sum(max(0,Population.cons),2)~=0;
      InfeasibleFrontNo = FrontNo(InfeasibleIndex);
      NextInfeasibleIndex = InfeasibleFrontNo<=MaxNo;
      InfeasibleIndex1 = InfeasibleIndex(NextInfeasibleIndex);
      Infeasible = Population(InfeasibleIndex1);
      InfeasibleFrontNo = FrontNo(InfeasibleIndex1); 
      Infeasible = InfeasibleUpdate(Infeasible,InfeasibleFrontNo,N);
      
end

function Population = FeasibleUpdate(Population,FrontNo,N)

if length(Population)>N
    % 给层排序，确定最大层
    % [FrontNo,MaxFNo] = NDSort(Population.objs,N);
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
    %最后获得的种群
    Population = [Population1,Population2]; 
end
end

function Population = InfeasibleUpdate(Population,FrontNo,N)

if length(Population) > N
    
    
    %% 给层排序，确定最大层
   % [FrontNo,MaxFNo] = NDSort([Population.objs,sum(max(0,Population.cons),2)],N);
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N);
    %% 确定前面的层
    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    
    %% 确定最后一层的个体和所需要的个体
    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    N1 = N- sum(Next); %所需要的个体
    %% 对最后一层进行Boudary sorting 
    [FrontNo1,MaxFNo1] = NDSort(-Population2.objs,N1);
   
    Next1 = FrontNo1 < MaxFNo1;
    Population3 = Population2(Next1);
    Last1 = FrontNo1==MaxFNo1;
    Population4 = Population(Last1);
   
    if sum(Last1)~= N1-sum(Next1)
    Population4 = Truncation(Population4,Population,N1-sum(Next1));
    end
    %% Population for next generation
    Population = [Population1,Population3,Population4];
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
