function Population = EnvironmentalSelection(Population,N)
% The environmental selection of NBA

%--------------------------------------------------------------------------
    %% Non-dominated sorting of the transformed problem
    [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);
    infeasible = sum(max(0,Population.cons),2)~=0;
    FrontNo = FrontNo(infeasible);
    Population = Population(infeasible);
    
  if size(Population,2)>N
     
      
      [FNSort,~]= sort(FrontNo);
      MaxFNo = FNSort(N); 
      Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front 
      Last     = find(FrontNo==MaxFNo);
      Population1= Population(Next);
    
    Population2= Population(Last);   
    N2 = N- sum(Next);
    %% Boudary sorting 
     [FrontNo2,MaxFNo2] = NDSort(-Population2.objs,N2);
     Next2 = FrontNo2 < MaxFNo2;
     
    %% Advanced crowding distance sorting
     CrowdDis = AdvancedCrowdingDistance(Population2.objs,sum(max(0,Population2.cons),2),FrontNo2);
     
    %% Select the solutions in the last front based on their crowding distances
     Last2     = find(FrontNo2==MaxFNo2);
    [~,Rank] = sort(CrowdDis(Last2),'descend');
    Next2(Last2(Rank(1:N2-sum(Next2)))) = true;
    
    %% Population for next generation
    Population2 = Population2(Next2);
    
    Population = [Population1,Population2];
    
  end
% 
%      if sum(Next)+size(Last,2)-N == 0
%          Next(Last)=1;
%      else
%      Del  = Truncation(Population,Last,sum(Next)+size(Last,2)-N);
%      Next(Last(~Del)) = true; 
%      end
%     % Temp = find(Last);
% 
%      
% 
%     %% Population for next generation
%     Population = Population(Next);
%    % FrontNo    = FrontNo(Next);
%    % CrowdDis   = CrowdDis(Next);
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
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(1./(1+CV(Front(Rank(j)))))*(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end

% function Population = EnvironmentalSelection(Population,N)
% % The environmental selection of NSGA-II
% 
% %--------------------------------------------------------------------------
% % Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% % research purposes. All publications which use this platform or any code
% % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% % for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% % Computational Intelligence Magazine, 2017, 12(4): 73-87".
% %--------------------------------------------------------------------------
% 
%     %% Non-dominated sorting
%     [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
%     Next = FrontNo < MaxFNo;
%     
%     %% Calculate the crowding distance of each solution
%    % CrowdDis = CrowdingDistance(Population.objs,FrontNo);
%     
%     %% Select the solutions in the last front based on their crowding distances
%     Last     = find(FrontNo==MaxFNo);
%     
%    % [~,Rank] = sort(CrowdDis(Last),'descend');
%    %Next(Last(Rank(1:N-sum(Next)))) = true;
%      if sum(Next)+size(Last,2)-N == 0
%          Next(Last)=1;
%      else
%      Del  = Truncation(Population,Last,sum(Next)+size(Last,2)-N);
%      Next(Last(~Del)) = true; 
%      end
%     % Temp = find(Last);
% 
%      
% 
%     %% Population for next generation
%     Population = Population(Next);
%    % FrontNo    = FrontNo(Next);
%    % CrowdDis   = CrowdDis(Next);
% end
% 
% 
% function Del = Truncation(Population,Last,K)
% % Select part of the solutions by truncation
% 
%     
% 
%     %% Truncation
%     Zmin       = min(Population.objs,[],1);
%    % Zmax = max(Population.objs,[],1);
%     PopObj = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(Population.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1)+1e-10)+1e-10;
%     PopObj = (PopObj(Last,:)); 
%     
%     Distance = pdist2(PopObj,PopObj);
%     Distance(logical(eye(length(Distance)))) = inf;
%     Del = false(1,size(PopObj,1));
%     while sum(Del) < K
%         Remain   = find(~Del);
%         Temp     = sort(Distance(Remain,Remain),2);
%         [~,Rank] = sortrows(Temp);
%         Del(Remain(Rank(1))) = true;
%     end
%       
% 
% end