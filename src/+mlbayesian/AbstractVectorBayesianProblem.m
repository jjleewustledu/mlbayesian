classdef AbstractVectorBayesianProblem < mlbayesian.IBayesianProblem 
	%% ABSTRACTVECTORBAYESIANPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a).  Copyright 2014 John Joowon Lee. 
 	%  $Id$ 
    
    properties         
        independentData % numeric
        dependentData   % numeric     
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
  		function this = AbstractVectorBayesianProblem(varargin) 
 			%% ABSTRACTVECTORBAYESIANPROBLEM 
 			%  Usage:  this = AbstractBayesianProblem([independent_data, dependent_data]) 
            
            p = inputParser;
            addOptional(p, 'indepData', [], @isnumeric);
            addOptional(p,   'depData', [], @isnumeric);
            parse(p, varargin{:});            
 			
            this.independentData           = p.Results.indepData;
            this.dependentData             = p.Results.depData;
            this.normSquaredDependentData_ = norm(this.dependentData)^2;
            assert(all(size(this.independentData) == size(this.dependentData)));
        end 
        function sse  = sumSquaredErrors(this, p)
            sse = norm(this.dependentData - this.estimateDataFast(p));
        end
        function x    = finalParams(this)
            x = this.bestFitParams;
        end 
    end
    
    properties (Access = 'protected')
        normSquaredDependentData_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

