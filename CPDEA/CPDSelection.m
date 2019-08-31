function [Population,Del,fCPD,LocalC] = CPDSelection(Population,K,eta,D_Dec,DominationX,Global)
% CPD based environmental selection in CPDEA

%--------------------------------------------------------------------------
% Copyright 2018-2019 Yiping Liu
% This is the code of CPDEA proposed in "Yiping Liu, Hisao Ishibuchi, Gary
% G. Yen, Yusuke Nojima and Naoki Masuyama, Handling Imbalance Between 
% Convergence and Diversity in the Decision Space in Evolutionary Multi-
% Modal Multi-Objective Optimization, IEEE Transactions on Evolutionary 
% Computation, 2019, Early Access, DOI: 10.1109/TEVC.2019.2938557".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Local Convergence Quality 
N = length(Population);
V = prod(Global.upper - Global.lower);
R = eta.*(V./N).^(1./Global.D);
temp = exp(-D_Dec.^2./(2.*R.^2))./(R.*(2.*pi).^.5); % Normal_distribution
LocalC = NaN(1,N);
for i = 1:N
    LocalC(i) = DominationX(i,:)*temp(:,i); 
end

%% Transform D_Dec based on Local Convergence Quality
D_t = D_Dec;
for i = 1:N-1
   for j = i+1:N
       x = (LocalC(i)+LocalC(j)).*0.5;       
       y = 1./((x+1).^1);            
       D_t(i,j) = D_t(i,j).*y; 
       D_t(j,i) = D_t(i,j);
   end
end
D_t(logical(eye(N))) = inf;

%% k-nearest neighbor
D_s = sort(D_t);
fCPD = 1./(1+sum(D_s(1:K,:)));

%% remove the worst solution
[~,Del] = max(fCPD);
Population(Del) = [];   
fCPD(Del) = [];
LocalC(Del) = [];   
end

   



