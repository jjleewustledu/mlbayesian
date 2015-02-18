classdef PolynomialProblem < mlbayesian.AbstractBayesianProblem
	%% POLYNOMIALPROBLEM is for testing IMCMC implementation

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
	methods 
        function this = PolynomialProblem(t, y)
            this = this@mlbayesian.AbstractBayesianProblem(t, y);
        end
        function this = estimateParameters(this)
            
            import mlbayesian.*;
            map = containers.Map;
            map('c0') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            map('c1') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            map('c2') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);       

            this.paramsManager = BayesianParameters(map);            
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end   
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('c0'), this.finalParams('c1'), this.finalParams('c2'));
        end
        function ed   = estimateDataFast(this, c0, c1, c2)
            ed = c0 + c1*this.independentData + c2*this.independentData.^2;
        end   
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

