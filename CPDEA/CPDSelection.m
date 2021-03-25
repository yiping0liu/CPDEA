function [Population,Del,fCPD,LocalC] = CPDSelection(Population,K,eta,D_Dec,DominationX,Prob)
% CPD based environmental selection in CPDEA

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Yen, G.G., Nojima, Y. and Masuyama, N., 2020. 
% Handling imbalance between convergence and diversity in the decision 
% space in evolutionary multimodal multiobjective optimization. IEEE 
% Transactions on Evolutionary Computation, 24(3), pp.551-565.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

%% Local Convergence Quality 
N = length(Population);
V = prod(Prob.upper - Prob.lower);
R = eta.*(V./N).^(1./Prob.D);
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

   



