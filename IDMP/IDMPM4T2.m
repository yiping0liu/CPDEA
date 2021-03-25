classdef IDMPM4T2 < PROBLEM
% <multi/many> <real> <multimodal>
% Imbalanced Distance Minimization Problems

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
        Points; % Vertexes
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M        = 4; 
            obj.D        = 4;
            obj.lower    = [-1,-1,-1,-1];
            obj.upper    = [1,1,1,1];
            obj.encoding = 'real';
            % Generate vertexes             
            pgon = nsidedpoly(obj.M);
            psize = [.1,.1,.1,.1];
            center = [-.5,-.5;.5,-.5;.5,.5;-.5,.5]; 
            obj.Points = NaN(obj.M,2,4);
            for i=1:4
                obj.Points(:,:,i) = pgon.Vertices.*psize(i)+center(i,:);
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N = size(PopDec,1);
            PopObj = NaN(N,obj.M);
            for i=1:obj.M
                temp = pdist2(PopDec(:,1:2),reshape(obj.Points(i,:,:),[2,4])');                
                temp(:,1) = temp(:,1)+100.*((abs(PopDec(:,3)+.6)).^2+(abs(PopDec(:,4)+.6)).^2);
                temp(:,2) = temp(:,2)+100.*((abs(PopDec(:,3)+.2)).^1.8+(abs(PopDec(:,4)+.2)).^1.8);
                temp(:,3) = temp(:,3)+100.*((abs(PopDec(:,3)-.2)).^1.6+(abs(PopDec(:,4)-.2)).^1.6);
                temp(:,4) = temp(:,4)+100.*((abs(PopDec(:,3)-.6)).^1.4+(abs(PopDec(:,4)-.6)).^1.4);
                PopObj(:,i) = min(temp,[],2);
            end             
        end
        
    end
end