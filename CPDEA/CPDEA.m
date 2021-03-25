classdef CPDEA < ALGORITHM
% <multi/many> <real> <multimodal>
% Handling Imbalance Between Convergence and Diversity in the Decision 
% Space in Evolutionary Multi-Modal Multi-Objective Optimization 
% K --- 3 --- parameter for k-nearest
% eta --- 2 --- parameter for local convergence

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Yen, G.G., Nojima, Y. and Masuyama, N., 2020. 
% Handling imbalance between convergence and diversity in the decision 
% space in evolutionary multimodal multiobjective optimization. IEEE 
% Transactions on Evolutionary Computation, 24(3), pp.551-565.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
   
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [K,eta] = Algorithm.ParameterSet(3,2);

            %% Initialization
            Population = Problem.Initialization(Problem.N+1); % randomly generate N+1 solutions    
            [Archive,fDKN] = ArchiveUpdate(Population,Problem.N,K); % update archive     
            D_Dec = pdist2(Population.decs,Population.decs,'euclidean'); % distance between pair of solutions in the decision space
            DominationX = zeros(Problem.N+1); % Pareto domination relationship between pair of solutions
            for i=1:Problem.N
               for j=i+1:Problem.N+1 
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
            [Population,Del,fCPD,~] = CPDSelection(Population,K,eta,D_Dec,DominationX,Problem); % Environmental Selection

            %% Optimization
            while Algorithm.NotTerminated(Archive)             
                %% Generation One Offspring
                p = 1; % probability of using the first reproduction operator
                if Problem.FE > Problem.maxFE*.5
                    p = 0.5;
                end
                sd = rand;
                if sd > p && length(Archive) > K % second reproduction operator
                    Parents = MatingSelection2(Archive,Population,fDKN,Problem.D);
                    Offspring = OperatorGAhalf(Parents);          
                else % first reproduction operator
                    MatingPool = TournamentSelection(2,2,fCPD);
                    Offspring = OperatorGAhalf(Population(MatingPool));
                end

                %% Update D_Dec and DominationX       
                D1 = pdist2(Offspring.decs,Population.decs,'euclidean');
                D1 = [D1(1:Del-1),0,D1(Del:end)];
                D_Dec(Del,:) = D1;
                D_Dec(:,Del) = D1';       
                OfpObj = repmat(Offspring.objs,Problem.N,1);
                L1 = OfpObj < Population.objs;
                L2 = OfpObj > Population.objs;
                LA = all((L1|(~L2))');
                LB = all((L2|(~L1))');
                LA = [LA(1:Del-1),0,LA(Del:end)];
                LB = [LB(1:Del-1),0,LB(Del:end)];
                DominationX(Del,:) = LB;
                DominationX(:,Del) = LA';

                %% Environmental Selection
                [Population,Del,fCPD,~] = CPDSelection([Population(1:Del-1),Offspring,Population(Del:end)],K,eta,D_Dec,DominationX,Problem); %%%%%%%%%%%

                %% Update Archive
                [Archive,fDKN] = ArchiveUpdate([Archive,Offspring],Problem.N,K); 
            end
        end
    end
end