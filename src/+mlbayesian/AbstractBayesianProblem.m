classdef (Abstract) AbstractBayesianProblem < mlbayesian.IBayesianProblem 
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
        expectedBestFitParams
        bestFitParams
        meanParams
        stdParams
    end
    
    methods %% GET
        function e = get.expectedBestFitParams(this)
            e = this.expectedBestFitParams_;
        end
        function p = get.bestFitParams(this)
            p = this.mcmc.bestFitParams;
        end
        function p = get.meanParams(this)
            p = this.mcmc.meanParams;
        end
        function p = get.stdParams(this)
            p = this.mcmc.stdParams;
        end
    end
    
    methods (Static)
        function idx  = indexOf(t, t0)
            %% INDEXOF finds the array-index closest to t0 for array t
            %  Usage:  idx = this.indexOf(t, t0)
            
            metric = abs(t - t0);
            [~,idx] = min(metric(:));
            idx = floor(idx) + 1;
        end
    end
    
	methods 
  		function this = AbstractBayesianProblem(varargin) 
 			%% ABSTRACTBAYESIANPROBLEM 
 			%  Usage:  this = AbstractBayesianProblem([independent_data, dependent_data]) 
            
            p = inputParser;
            addOptional(p, 'indepData', [], @isnumeric); %% && this.uniformSampling(x)
            addOptional(p,   'depData', [], @isnumeric);
            parse(p, varargin{:});            
 			
            this.independentData = this.offsetZeros(p.Results.indepData);
            this.dependentData   = p.Results.depData;
            assert(all(size(this.independentData) == size(this.dependentData)));
        end 
        function sse  = sumSquaredErrors(this, p)
            p   = num2cell(p);        
            sse = sum(abs(this.dependentData - this.estimateDataFast(p{:})).^2);
            if (sse < eps)
                sse = sse + (1 + rand(1))*eps; 
            end
            %assert(isfinite(sse) && ~isnan(sse), 'AbstractBayesianProblem.p -> %s', cell2str(p));
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
        function x    = finalMeans(this, key)
            x = this.meanParams(this.paramsManager.paramsIndices(key));
        end   
        function x    = finalStds(this, key)
            x = this.stdParams(this.paramsManager.paramsIndices(key));
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
    
    %% PROTECTED
    
    properties (Access = 'protected')
        expectedBestFitParams_
    end
    
    methods (Access = 'protected')
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
            tf  = all(abs(dt1*ones(1,length(dts)) - dts) < eps('single'));
        end
        function t  = ensureRow(~, t)
            if (size(t,1) > size(t,2))
                t = t'; 
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

