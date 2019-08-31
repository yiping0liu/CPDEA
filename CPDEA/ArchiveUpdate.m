function [Population,fDKN] = ArchiveUpdate(Population,N,K)
% Archive Update in CPDEA

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

    FrontNo = NDSort(Population.objs,inf);   
    Population = Population(FrontNo ==1);
    
    [Choose,fDKN] = DoubleNearestSelection(Population.objs,Population.decs,N,K);
    
    Population = Population(Choose);
    fDKN = fDKN(Choose);

end

function [Choose,fDN] = DoubleNearestSelection(PopObj,PopDec,N,K)
% Select solutions based on double k-nearest neighbor method

    Np = size(PopObj,1);
    Choose = true(1,Np);
    
    if Np <= K
        fDN = zeros(1,Np);
        return;
    end   
    
    d_obj = pdist2(PopObj,PopObj,'euclidean');
    d_dec = pdist2(PopDec,PopDec,'euclidean');
    d_obj(logical(eye(Np))) = inf;
    d_dec(logical(eye(Np))) = inf;
    
    sdo = sort(d_obj);
    sdd = sort(d_dec);    
    dn_obj = sum(sdo(1:K,:));
    dn_dec = sum(sdd(1:K,:));
    avg_dn_obj = mean(dn_obj);
    avg_dn_dec = mean(dn_dec);
    if avg_dn_obj == 0
        avg_dn_obj = inf;
    end
    if avg_dn_dec == 0
        avg_dn_dec = inf;
    end
    fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
    
                
    while sum(Choose) > N
        [~,Del] = max(fDN);
        Choose(Del) = false;
        d_obj(Del,:) = inf;
        d_obj(:,Del) = inf;
        d_obj(Del,:) = inf;
        d_obj(:,Del) = inf;
        
        sdo = sort(d_obj);
        sdd = sort(d_dec);    
        dn_obj = sum(sdo(1:K,:));
        dn_dec = sum(sdd(1:K,:));
        avg_dn_obj = mean(dn_obj);
        avg_dn_dec = mean(dn_dec);
        if avg_dn_obj == 0
            avg_dn_obj = inf;
        end
        if avg_dn_dec == 0
            avg_dn_dec = inf;
        end
        fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
               
        fDN(~Choose) = -inf;
    end         
end