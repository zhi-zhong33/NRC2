 function Parents =MatingSelection(Feasible,Infeasible,N)
 
 
 
    Parents = [];
    %% 确定比较指标 可行解根据多样性，不可行解根据
    
    Population = [Feasible,Infeasible];
    [FrontNo,~] = NDSort([Population.objs,sum(max(0,Population.cons),2)],inf);
    Density = DensityEstimate(Population);  
%     [~,rank]= sort(sum(max(0,Population.cons),2));
%     CV = rank./max(rank);
%     FitAll = FrontNo + CV + Density;
     FitAll =  Density;
   
    %% 根据概率进行选择
        IndexInf1 = floor(rand(1,N)*length(Population))+1;
        IndexInf2 = floor(rand(1,N)*length(Population))+1; 

         for i = 1:N  
            if FitAll(IndexInf1(i))<= FitAll(IndexInf2(i))
               ParentTemp = Population(IndexInf1(i));
           else
               ParentTemp = Population(IndexInf2(i));
            end
            Parents = [Parents,ParentTemp];
         end  

 end
 
 function Density = DensityEstimate(Population)
    PopObjNew = Population.objs;
  % PopObjNew = Population.decs;
    Zmin       = min(PopObjNew,[],1)-1e-10;
    Zmax = max(PopObjNew,[],1);
    PopObj = (PopObjNew-repmat(Zmin,length(Population.objs),1))./(repmat(Zmax,length(Population.objs),1)-repmat(Zmin,length(Population.objs),1));
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    DistanceSort = sort(Distance,2);
    DensityValue = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));
    [~,rank]= sort(DensityValue);
%     DensityCV =  sum(max(0,Population.cons),2)./max(sum(max(0,Population.cons),2));
%     DensityValue(sum(max(0,Population.cons),2)~=0) = DensityCV(sum(max(0,Population.cons),2)~=0);
   % DensityValue    = min(Distance,[],2);
    Density =  rank./max(rank);
end