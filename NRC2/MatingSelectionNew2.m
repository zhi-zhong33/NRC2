 function Parents =MatingSelectionNew2(Feasible,Infeasible,N)
    %% 确定比较指标 可行解根据多样性，不可行解根据
    Parents = [];
    Population = [Feasible,Infeasible];
   [FrontNo1,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);
  % [FrontNo1,~] = NDSort(Population.objs,inf);
  %  BoundaryNo = BSort(Population,FrontNo1,length(Population));
  % [FrontNo2,~] = NDSort(Population.objs,sum(max(0,Population.cons),2),inf);
    DensityValue = DensityCal(Population);%[DensityCal(Feasible);DensityCal(Infeasible)];
     Fitness1= FrontNo1+DensityValue;
  %  Fitness1 = FrontNo1+0.1*BoundaryNo+0.01*DensityValue;
   % Fitness1 = DensityValue;
   % Fitness2 = FrontNo1+0.1*BoundaryNo+0.01*DensityValue;
     Fitness2 = FrontNo1+DensityValue;
    
        Index1 = floor(rand(1,N)*length(Population))+1;
        Index2 = floor(rand(1,N)*length(Population))+1; 
         for i = 1:2:N    
                     if Fitness1(Index1(i))< Fitness1(Index2(i))
                         ParentTemp1 = Population(Index1(i));
                     else
                         ParentTemp1 = Population(Index2(i));
                     end         
               
                    if Fitness2(Index1(i+1))< Fitness2(Index2(i+1))
                        ParentTemp2 = Population(Index1(i+1));
                    else
                        ParentTemp2 = Population(Index2(i+1));
                    end                 
            Parents = [Parents,ParentTemp1,ParentTemp2];
         end  

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
 
 
 
 function DensityValue =DensityEstimatNew(Population,FrontNo,BoundaryNo,CV)
 Index= CV==0;
 Feasible = Population(Index);
 FrontNoF = FrontNo(Index);
 BoundaryNoF=BoundaryNo(Index);
 DensityValue =  zeros(1,length(Population));
 if ~isempty(Feasible)
 DensityTempF = DensityEstimate(Feasible,FrontNoF,BoundaryNoF,length(Feasible));
 DensityValue(Index)=DensityTempF;
 end
 
 Index2= CV==0;
 Infeasible = Population(Index2);
 FrontNoInf = FrontNo(Index2);
 BoundaryNoInf=BoundaryNo(Index2);
 if ~isempty(Infeasible)
 DensityTempF = DensityEstimate(Infeasible,FrontNoInf,BoundaryNoInf,length(Infeasible));
 DensityValue(Index2)=DensityTempF;
 end
 end
 
 
 function DensityValue = DensityEstimate(Population,FrontNo,BoundaryNo,N)
  DensityValue = zeros(N,1);
  Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        BoundaryTemp = BoundaryNo(Front);
        PopulationTemp = Population(Front);
        DensityValueTemp = DensityValue(Front);
        Boundarys = setdiff(unique(BoundaryTemp),inf);
        for i = 1:length(Boundarys)
            Boundary =  find(BoundaryTemp==Boundarys(i));
            DensityTemp = DensityCal(PopulationTemp(Boundary)); 
           % DensityTemp = CrowdingDistance(PopulationTemp(Boundary).objs,Boundary);
            DensityValueTemp(Boundary)=DensityTemp;
        end 
        DensityValue(Front) = DensityValueTemp;   
    end
 end
 
 
%  function CrowdDis = CrowdingDistance(PopObj,FrontNo)
% % Calculate the crowding distance of each solution front by front
%     [N,M]    = size(PopObj);
%     CrowdDis = zeros(1,N);
%     Fronts   = setdiff(unique(FrontNo),inf);
%     for f = 1 : length(Fronts)
%         Front = find(FrontNo==Fronts(f));
%         Fmax  = max(PopObj(Front,:),[],1);
%         Fmin  = min(PopObj(Front,:),[],1);
%         for i = 1 : M
%             [~,Rank] = sortrows(PopObj(Front,i));
%             CrowdDis(Front(Rank(1)))   = inf;
%             CrowdDis(Front(Rank(end))) = inf;
%             for j = 2 : length(Front)-1
%               % CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(1./(1+CV(Front(Rank(j)))))*(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
%               CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
%             end
%         end
%     end
% end
 
 
 function Density = DensityCal(Population)
 
    if isempty(Population)
        Density =[];
    elseif length(Population) <= 2
        Density = zeros(length(Population),1);
    else
    Zmin       = min(Population.objs,[],1);
    PopObj = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(Population.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1));
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = Inf;
    DistanceSort = sort(Distance,2);
    Density = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));
 %   Density = 1./(1+DistanceSort(:,floor(sqrt(N))+2));
 %   [~,rank]= sort(DensityValue);
%     DensityCV =  sum(max(0,Population.cons),2)./max(sum(max(0,Population.cons),2));
%     DensityValue(sum(max(0,Population.cons),2)~=0) = DensityCV(sum(max(0,Population.cons),2)~=0);
  %  Density    = 1./(1+min(Distance,[],2));
 %   Density =  rank./max(rank);
    end
 end
% 
%   function Density = DensityCal(Population,N)
%  
%     if isempty(Population)
%         Density =[];
%     elseif length(Population) <= 2
%         Density = zeros(length(Population),1);
%     else
%         PopObjTemp = [Population.objs,Population.cons];     
%         Zmin       = min(PopObjTemp,[],1);
%         PopObj = (PopObjTemp-repmat(Zmin,length(Population.objs),1))./(repmat(max(PopObjTemp),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1)+1e-10)+1e-10;
%        Distance = pdist2(PopObj,PopObj);
%        Distance(logical(eye(length(Distance)))) = inf;
%        DistanceSort = sort(Distance,2);
%     %Density = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));
%     Density = 1./(1+DistanceSort(:,floor(sqrt(N))+1));
%  %   [~,rank]= sort(DensityValue);
% %     DensityCV =  sum(max(0,Population.cons),2)./max(sum(max(0,Population.cons),2));
% %     DensityValue(sum(max(0,Population.cons),2)~=0) = DensityCV(sum(max(0,Population.cons),2)~=0);
%   %  Density    = 1./(1+min(Distance,[],2));
%  %   Density =  rank./max(rank);
%     end
% end