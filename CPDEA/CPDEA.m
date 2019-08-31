function CPDEA(Global)
% <algorithm> <A-G>
% Handling Imbalance Between Convergence and Diversity in the Decision 
% Space in Evolutionary Multi-Modal Multi-Objective Optimization 
% K --- 3 --- parameter for k-nearest
% eta --- 2 --- parameter for local convergence

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

    %% Parameter setting
    [K,eta] = Global.ParameterSet(3,2);
     
    %% Initialization
    Population = Global.Initialization(Global.N+1); % randomly generate N+1 solutions    
    [Archive,fDKN] = ArchiveUpdate(Population,Global.N,K); % update archive     
    D_Dec = pdist2(Population.decs,Population.decs,'euclidean'); % distance between pair of solutions in the decision space
    DominationX = zeros(Global.N+1); % Pareto domination relationship between pair of solutions
    for i=1:Global.N
       for j=i+1:Global.N+1 
           L1 = Population(i).objs < Population(j).objs;
           L2 = Population(i).objs > Population(j).objs;
           if all(L1|(~L2))
               DominationX(i,j) = 0;
               DominationX(j,i) = 1;
           elseif all(L2|(~L1))
               DominationX(i,j) = 1;
               DominationX(j,i) = 0;
           end           
       end
    end   
    [Population,Del,fCPD,~] = CPDSelection(Population,K,eta,D_Dec,DominationX,Global); % Environmental Selection
      
    %% Optimization
    while Global.NotTermination(Archive)             
        %% Generation One Offspring
        p = 1; % probability of using the first reproduction operator
        if Global.evaluated > Global.evaluation*.5
            p = 0.5;
        end
        sd = rand;
        if sd > p && length(Archive) > K % second reproduction operator
            Parents = MatingSelection2(Archive,Population,fDKN,Global);
            Offspring = Global.Variation(Parents,1);            
        else % first reproduction operator
            MatingPool = TournamentSelection(2,2,fCPD);
            Offspring = Global.Variation(Population(MatingPool),1);
        end
                       
        %% Update D_Dec and DominationX       
        D1 = pdist2(Offspring.decs,Population.decs,'euclidean');
        D1 = [D1(1:Del-1),0,D1(Del:end)];
        D_Dec(Del,:) = D1;
        D_Dec(:,Del) = D1';       
        OfpObj = repmat(Offspring.objs,Global.N,1);
        L1 = OfpObj < Population.objs;
        L2 = OfpObj > Population.objs;
        LA = all((L1|(~L2))');
        LB = all((L2|(~L1))');
        LA = [LA(1:Del-1),0,LA(Del:end)];
        LB = [LB(1:Del-1),0,LB(Del:end)];
        DominationX(Del,:) = LB;
        DominationX(:,Del) = LA';
                
        %% Environmental Selection
        [Population,Del,fCPD,~] = CPDSelection([Population(1:Del-1),Offspring,Population(Del:end)],K,eta,D_Dec,DominationX,Global);
        
        %% Update Archive
        [Archive,fDKN] = ArchiveUpdate([Archive,Offspring],Global.N,K);                
              
    end
end