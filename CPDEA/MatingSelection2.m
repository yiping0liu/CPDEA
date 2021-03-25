function Parents = MatingSelection2(Archive,Population,fDKN,D)
% Mating selection for the second reproduction operator in CPDEA

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Yen, G.G., Nojima, Y. and Masuyama, N., 2020. 
% Handling imbalance between convergence and diversity in the decision 
% space in evolutionary multimodal multiobjective optimization. IEEE 
% Transactions on Evolutionary Computation, 24(3), pp.551-565.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

% Parent 1
mfDKN = min(fDKN);
CP1 = find(fDKN == mfDKN);
P1 = CP1(randi(length(CP1)));
Parent1 = Archive(P1);

% Parent 2
N2 = D;
Pop2 = [Archive,Population];
d_dec = pdist2(Parent1.decs,Pop2.decs,'euclidean');
[~,I] = sort(d_dec);
I = I(2:N2+1);
P2 = I(randi(N2));

Parents = [Parent1,Pop2(P2)];
end