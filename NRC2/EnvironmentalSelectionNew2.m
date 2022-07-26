function [Feasible,Infeasible] = EnvironmentalSelectionNew2(Population,N)
      
      % ȥ���ظ���
      [~,b]=unique(Population.objs,'rows');
      if size(b,1)> N
           Population = Population(1,b);
      end
      % ��Լ������һ��Ŀ�꣬�ֲ�
      [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);
      
      % ��ȡ���еĿ��н⣬���һ�ÿ��н�ı��
      FeasibleIndex = sum(max(0,Population.cons),2)==0;
      Feasible = Population(FeasibleIndex);
      FrontNoFeasible= FrontNo(FeasibleIndex);
      Feasible = FeasibleUpdate(Feasible,FrontNoFeasible,N);
      
      % ��ȡ�����н⣬���һ����Ӧ���
      InfeasibleIndex = sum(max(0,Population.cons),2)~=0;
      Infeasible = Population(InfeasibleIndex);
      FrontNoInfeasible= FrontNo(InfeasibleIndex);
      Infeasible = InfeasibleUpdate(Infeasible,FrontNoInfeasible,N);
      
      
   %   Ratio =  sum(FrontNoFeasible==1)./(sum(FrontNoFeasible==1)+sum(FrontNoInfeasible==1));  
end

function Population = FeasibleUpdate(Population,FrontNo,N)

if length(Population)>N
    % ��������ȷ������
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
    %ȷ�����Ĳ�
    Population = [Population1,Population2]; 
end
end

function Population = InfeasibleUpdate(Population,FrontNo,N)

if length(Population) > N
    %% ��������ȷ������
    FrontNoSort = sort(FrontNo);
	BoundaryNo = BSort(Population,FrontNo,length(Population));
	DensityValue =  DensityCal(Population);%[DensityCal(Feasible);DensityCal(Infeasible)];
    Fitness = FrontNo+0.1*BoundaryNo./max(BoundaryNo)+0.01*DensityValue./max(DensityValue);
    [~,index]=sort(Fitness);
	Population = Population(index(1:N));
    
    
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


 function  BoundaryNo = BSort(Population,FrontNo,N)
 BoundaryNo = zeros(1,N);
 Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        [BoundaryNoTemp,~]= NDSort(-Population(Front).objs,inf);
        BoundaryNo(Front)=BoundaryNoTemp;
    end
 end


 function Density = DensityCal(Population)
    if isempty(Population)
        Density=[];
    elseif length(Population)<=2
        Density = zeros(length(Population),1);
    else
    Zmin       = min(Population.objs,[],1);
    PopObj = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(Population.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1)+1e-10)+1e-10;
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    DistanceSort = sort(Distance,2);
    Density = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));
 %   [~,rank]= sort(DensityValue);
%     DensityCV =  sum(max(0,Population.cons),2)./max(sum(max(0,Population.cons),2));
%     DensityValue(sum(max(0,Population.cons),2)~=0) = DensityCV(sum(max(0,Population.cons),2)~=0);
   %  Density    = 1./(1+min(Distance,[],2));
 %   Density =  rank./max(rank);
    end
 end
