function Parents = MatingSelection2(Archive,Population,fDKN,Global)
% Mating selection for the second reproduction operator in CPDEA

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

% Parent 1
mfDKN = min(fDKN);
CP1 = find(fDKN == mfDKN);
P1 = CP1(randi(length(CP1)));
Parent1 = Archive(P1);

% Parent 2
N2 = Global.D;
Pop2 = [Archive,Population];
d_dec = pdist2(Parent1.decs,Pop2.decs,'euclidean');
[~,I] = sort(d_dec);
I = I(2:N2+1);
P2 = I(randi(N2));

Parents = [Parent1,Pop2(P2)];
end