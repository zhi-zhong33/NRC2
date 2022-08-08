function [Feasible,Infeasible] = EnvironmentalSelectionElite(Population,N)
      
      % remove duplicated solutions
      [~,b]=unique(Population.objs,'rows');
      Population = Population(1,b);
      
      [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);

      % find the feasible soltuons and the corresponding indexes
      FeasibleIndex = sum(max(0,Population.cons),2)==0;    
      Feasible = Population(FeasibleIndex);
      FeasibleFrontNo = FrontNo(FeasibleIndex);
      Feasible = FeasibleUpdate(Feasible,FeasibleFrontNo,N); 
% 
      if length(Feasible)<N 
         N1 = N-length(Feasible);
         InfeasibleIndex = sum(max(0,Population.cons),2)~=0;
         Infeasible = Population(InfeasibleIndex); 
         [~,Index]=sort(sum(max(0,Infeasible.cons),2));
         Feasible=[Feasible,Infeasible(Index(1:N1))];
      end
      
      
      % find the infeasible solutions and the coresponding indexes
      
         InfeasibleIndex = sum(max(0,Population.cons),2)~=0;
         Infeasible = Population(InfeasibleIndex);
         InfeasibleFrontNo = FrontNo(InfeasibleIndex);
         Temp = InfeasibleFrontNo==1;
         Infeasible =Infeasible(Temp);
         InfeasibleFrontNo = InfeasibleFrontNo(Temp);
         Infeasible = InfeasibleUpdate(Infeasible,InfeasibleFrontNo, N );
 
end

function Population = FeasibleUpdate(Population,FrontNo,N)

if length(Population)>N
  
    % [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N);
    

    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    

    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    if sum(Last)~= N-sum(Next)
    Population2 = Truncation(Population2,Population,N-sum(Next));
    end
  
    Population = [Population1,Population2]; 
end
end

function Population = InfeasibleUpdate(Population,FrontNo,N)

if length(Population) > N
    
    

 
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N);

    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    

    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    N1 = N- sum(Next); 
    %% Boudary sorting 
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
    PopObjTemp = Population.objs;
    PopObj = (PopObjTemp -repmat(Zmin,length(Population),1))./(repmat(max(PopAll.objs),length(Population),1)-repmat(Zmin,length(Population),1));
 
    
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
