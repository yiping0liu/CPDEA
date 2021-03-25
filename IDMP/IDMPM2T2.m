classdef IDMPM2T2 < PROBLEM
% <multi/many> <real> <multimodal>
% Imbalanced Distance Minimization Problems
% a --- 0.4 --- alpha

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Yen, G.G., Nojima, Y. and Masuyama, N., 2020. 
% Handling imbalance between convergence and diversity in the decision 
% space in evolutionary multimodal multiobjective optimization. IEEE 
% Transactions on Evolutionary Computation, 24(3), pp.551-565.
%------------------------------- Copyright --------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    properties(Access = private)
        a;      % alpha
        Points; % Vertexes
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M        = 2; 
            obj.D        = 2;
            obj.a       = obj.ParameterSet(.4);
            obj.lower    = [-1,-1];
            obj.upper    = [1,1];
            obj.encoding = 'real';
            % Generate vertexes             
            psize = 0.10.*ones(1,2);
            center = [-0.50,0.50];
            obj.Points = [center - psize; center + psize];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N = size(PopDec,1);
            PopObj = NaN(N,obj.M);
            for i=1:obj.M               
                temp = abs(repmat(PopDec(:,1),1,2)-repmat(obj.Points(i,:),N,1));
                temp(:,1) = temp(:,1)+100.*(abs(PopDec(:,2)+0.5)).^2;
                temp(:,2) = temp(:,2)+100.*(abs(PopDec(:,2)-0.5)).^(2-obj.a);                
                PopObj(:,i) = min(temp,[],2);
            end   
        end
        
    end
end