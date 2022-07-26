 function Parents =MatingSelectionNew(Feasible,Infeasible,N)
    %% 确定比较指标 可行解根据多样性，不可行解根据
    Parents = [];
    Population = [Feasible,Infeasible];
    [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);
    BoundaryNo = BSort(Population,FrontNo,length(Population));
    CV = sum(max(0,Population.cons),2);
   % DensityValue = DensityEstimatNew(Population,FrontNo,BoundaryNo,CV);
   DensityValue = DensityEstimatNew(Population,FrontNo,BoundaryNo,CV);
   % DensityValue =DensityEstimate(Population,FrontNo,BoundaryNo,length(Population)); 
    %DensityValue = DensityValue./(1+CV'); 
%    Index= CV~=0;
%    CV(Index)=1;
   %MatingPool = TournamentSelection(2,length(Population),CV,FrontNo,BoundaryNo,DensityValue);
  % Parents = Population(MatingPool(1:N));
 %  FitAll = FrontNo+0.1*BoundaryNo;
        Index1 = floor(rand(1,N)*length(Population))+1;
        Index2 = floor(rand(1,N)*length(Population))+1; 
         for i = 1:N  
            if FrontNo(Index1(i))< FrontNo(Index2(i))
               ParentTemp = Population(Index1(i));
            elseif FrontNo(Index1(i))> FrontNo(Index2(i))
               ParentTemp = Population(Index2(i));
            else
                if rand < 0.5
                     if BoundaryNo(Index1(i))< BoundaryNo(Index2(i)) || ( BoundaryNo(Index1(i))== BoundaryNo(Index2(i))&& CV(Index1(i))< CV(Index2(i))) || ( BoundaryNo(Index1(i))== BoundaryNo(Index2(i))&& CV(Index1(i))==CV(Index2(i)) &&  DensityValue(Index1(i))< DensityValue(Index2(i)))
                        ParentTemp = Population(Index1(i));
                    else
                        ParentTemp = Population(Index2(i));
                     end         
                else
                    if DensityValue(Index1(i))< DensityValue(Index2(i))
                        ParentTemp = Population(Index1(i));
                    else
                        ParentTemp = Population(Index2(i));
                    end              
                end
                    
                    
%                 if rand <0.5
%                     if CV(Index1(i))< CV(Index2(i)) || ( CV(Index1(i))== CV(Index2(i))&& DensityValue(Index1(i))< DensityValue(Index2(i)))
%                         ParentTemp = Population(Index1(i));
%                     else
%                         ParentTemp = Population(Index2(i));
%                     end
%                 else
%                     if DensityValue(Index1(i))< DensityValue(Index2(i))
%                         ParentTemp = Population(Index1(i));
%                     else
%                         ParentTemp = Population(Index2(i));
%                     end
%                 end
                    
            end
            Parents = [Parents,ParentTemp];
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
 
   
    if length(Population) <= 2
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