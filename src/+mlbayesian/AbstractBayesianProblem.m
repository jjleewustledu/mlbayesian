classdef AbstractBayesianProblem < mlbayesian.IBayesianProblem 
	%% ABSTRACTBAYESIANPROBLEM  
    %  Yet abstract:   
    %      properties showPlots
    %      methods estimateParameters, estimateData, estimateDataFast

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties         
        independentData % numeric, e.g., times
        dependentData   % numeric, e.g., densities = f(time)
        paramsManager   % mlbayesian.IBayesianParameters
        mcmc            % mlbayesian.IMCMC
    end
    
    properties (Dependent)
        bestFitParams
    end
    
    methods %% GET
        function p = get.bestFitParams(this)
            p = this.mcmc.bestFitParams;
        end
    end
    
	methods 
  		function this = AbstractBayesianProblem(varargin) 
 			%% ABSTRACTBAYESIANPROBLEM 
 			%  Usage:  this = AbstractBayesianProblem([independent_data, dependent_data]) 
            
            p = inputParser;
            addOptional(p, 'indepData', [], @isnumeric);
            addOptional(p,   'depData', [], @isnumeric);
            parse(p, varargin{:});            
 			
            this.independentData = p.Results.indepData;
            this.dependentData   = p.Results.depData;
            assert(all(size(this.independentData) == size(this.dependentData)));
        end 
        function sse  = sumSquaredErrors(this, p)
            p   = num2cell(p);
            sse = norm(this.dependentData - this.estimateDataFast(p{:}));
        end
        function x    = finalParams(this, key)
            x = this.bestFitParams(this.paramsManager.paramsIndices(key));
        end 
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

