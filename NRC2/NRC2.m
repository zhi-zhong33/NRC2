classdef NRC2 < ALGORITHM
% <multi> <real> <constrained>
% NRC2
% 
%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
 methods
    function main(Algorithm,Problem)
    Population = Problem.Initialization();
    ArcPop = []; 
    %% Optimization
    while Algorithm.NotTerminated(Population)
             Parents = MatingSelectionElite(Population,ArcPop,Problem.N);
             Offspring  = OperatorGA(Parents); 
             [Population,ArcPop] = EnvironmentalSelectionElite([Population,ArcPop,Offspring],Problem.N);      
    end
    end
 end
end