function [Feasible,Infeasible] = EnvironmentalSelectionPro(Population,N)
      
      % ȥ���ظ���
      [~,b]=unique(Population.objs,'rows');
      Population = Population(1,b);
      
      [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);

      % ��ȡ���еĿ��н⣬���һ�ÿ��н�ı��
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
      
      % ��ȡ�����н⣬���һ����Ӧ��� 
    
         InfeasibleIndex = sum(max(0,Population.cons),2)~=0;
         Infeasible = Population(InfeasibleIndex);
         InfeasibleFrontNo = FrontNo(InfeasibleIndex);
         Infeasible = InfeasibleUpdate(Infeasible,InfeasibleFrontNo, N ); 
 
end

function Population = FeasibleUpdate(Population,FrontNo,N)

if length(Population)>N
    % ��������ȷ������
    % [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N);
    
    % ȷ��ǰ��Ĳ�
    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    
    % ȷ�����һ�㲢ɾ������ĸ���
    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    if sum(Last)~= N-sum(Next)
    Population2 = Truncation(Population2,Population,N-sum(Next));
    end
    %����õ���Ⱥ
    Population = [Population1,Population2]; 
end
end

function Population = InfeasibleUpdate(Population,FrontNo,N)

if length(Population) > N
    
    
    %% ��������ȷ������
   % [FrontNo,MaxFNo] = NDSort([Population.objs,sum(max(0,Population.cons),2)],N);
    FrontNoSort = sort(FrontNo);
    MaxFNo = FrontNoSort(N);
    %% ȷ��ǰ��Ĳ�
    Next = FrontNo < MaxFNo;
    Population1 = Population(Next);
    
    %% ȷ�����һ��ĸ��������Ҫ�ĸ���
    Last = FrontNo == MaxFNo;
    Population2 = Population(Last);
    N1 = N- sum(Next); %����Ҫ�ĸ���
    %% �����һ�����Boudary sorting 
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
