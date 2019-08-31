function varargout = IDMPM2T2(Operation,Global,input)
% <problem> <IDMP>
% Imbalanced Distance Minimization Problems
% a --- 0.4 --- alpha
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright 2018-2019 Yiping Liu
% This is the code of benchmarks proposed in "Yiping Liu, Hisao Ishibuchi, 
% Gary G. Yen, Yusuke Nojima and Naoki Masuyama, Handling Imbalance Between 
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
    a = Global.ParameterSet(.4);
    persistent Points;
    switch Operation
        case 'init'
            Global.M          = 2;
            Global.D          = 2;
            Global.M          = 2;
            Global.D          = 2;
            Global.lower      = [-1,-1];
            Global.upper      = [1,1];
            Global.operator   = @EAreal;           
            NP = 2;
            psize = 0.10.*ones(1,NP);
            center = [-0.50,0.50];
            Points = [center - psize; center + psize];          
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            N = size(PopDec,1);
            PopObj = NaN(N,Global.M);           
            NP = 2;                      
            for i=1:Global.M
                temp = abs(repmat(PopDec(:,1),1,NP)-repmat(Points(i,:),N,1));
                temp(:,1) = temp(:,1)+100.*(abs(PopDec(:,2)+0.5)).^2;
                temp(:,2) = temp(:,2)+100.*(abs(PopDec(:,2)-0.5)).^(2-a);                
                PopObj(:,i) = min(temp,[],2);
            end                        
            PopCon = [];            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            div = 1/input;
            temp = 0:div:0.2;
            f = [temp',0.2-temp'];
            varargout = {f};
    end
end