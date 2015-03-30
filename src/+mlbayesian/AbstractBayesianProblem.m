classdef AbstractBayesianProblem < mlbayesian.IBayesianProblem 
	%% ABSTRACTBAYESIANPROBLEM  
    %  Yet abstract:   
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
    
    methods (Static)        
        function idx  = indexOf(t, t0)
            metric = abs(t - t0);
            [~,idx] = min(metric(:));
            idx = idx + 1;
        end
    end
    
	methods 
  		function this = AbstractBayesianProblem(varargin) 
 			%% ABSTRACTBAYESIANPROBLEM 
 			%  Usage:  this = AbstractBayesianProblem([independent_data, dependent_data]) 
            
            p = inputParser;
            addOptional(p, 'indepData', [], @(x) isnumeric(x) && this.uniformSampling(x));
            addOptional(p,   'depData', [], @isnumeric);
            parse(p, varargin{:});            
 			
            this.independentData = this.offsetZeros(p.Results.indepData);
            this.dependentData   = p.Results.depData;
            assert(all(size(this.independentData) == size(this.dependentData)));
        end 
        function sse  = sumSquaredErrors(this, p)
            p   = num2cell(p);
            sse = sum(abs(this.dependentData - this.estimateDataFast(p{:})).^2);
        end
        function ps   = adjustParams(~, ps)
            %% ADJUSTPARAMS:  override as needed for parameter constraints
        end
        function q    = Q(this)
            q = this.sumSquaredErrors(this.bestFitParams);
        end
        function nq   = normalizedQ(this)
            nq = this.Q/sum(abs(this.dependentData).^2);
        end
        function x    = finalParams(this, key)
            x = this.bestFitParams(this.paramsManager.paramsIndices(key));
        end      
        function ensureKeyOrdering(this, currentKeys)
            storedKeys = this.paramsManager.paramsMap.keys;
            for k = 1:length(storedKeys)
                assert(strcmp(storedKeys{k}, currentKeys{k}), ...
                       sprintf('AbstractBayesianProblem.ensureKeyOrdering:  expected %s but received %s', ...
                       storedKeys{k}, currentKeys{k}));
            end
        end
    end
    
    %% PRIVATE
    
    methods (Access = 'private')
        function x = ditherZeros(~, x)
            if (any(0 == x))
                x = x + (rand(1) + 1) * eps;
            end
        end
        function x = offsetZeros(~, x)            
            if (any(0 == x))
                x = x + eps;
            end
        end
        function tf = uniformSampling(this, t)
            t   = this.ensureRow(t);
            dts = t(2:end) - t(1:end-1);
            dt1 = t(2) - t(1);
            tf  = all(dt1*ones(1,length(dts)) == dts);
        end
        function t  = ensureRow(~, t)
            if (size(t,1) > size(t,2))
                t = t'; 
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

