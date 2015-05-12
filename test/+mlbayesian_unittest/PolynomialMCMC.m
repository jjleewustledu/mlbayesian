classdef PolynomialMCMC < mlbayesian.AbstractMcmcProblem
	%% POLYNOMIALMCMC is for testing MCMC implementations; MCMC requires access to an AbstractMcmcProblem

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties 
        baseTitle = 'no title'
        xLabel    = 'x label'
        yLabel    = 'y label'
    end
    
%     properties
%         independentData
%         dependentData
%         paramsManager
%         mcmc    
%     end
%     
%     properties (Dependent)
%         bestFitParams
%     end
%     
%     methods %% GET
%         function p = get.bestFitParams(this)
%             p = this.mcmc.bestFitParams;
%         end
%     end
    
	methods
        function this = PolynomialMCMC(indDat, depDat)
            this.independentData = indDat;
            this.dependentData = depDat;
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
        function sse  = sumSquaredErrors(this, p)
            p   = num2cell(p);
            sse = norm(this.dependentData - this.estimateDataFast(p{:}));
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.bestFitParams(1), this.bestFitParams(2), this.bestFitParams(3));
        end
        function ed   = estimateDataFast(this, c0, c1, c2)
            ed = c0 + c1*this.independentData + c2*this.independentData.^2;
        end   
        function x    = finalParams(this, key)
            x = this.bestFitParams(this.paramsManager.paramsIndices(key));
        end 
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

