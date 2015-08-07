classdef BayesianParameters < mlbayesian.IBayesianParameters
	%% BAYESIANPARAMETERS 
    
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$  	 
    
    properties        
        nProposals = 100   % number of loops in parameter prob phase
        nPop       =  50   % number of population
        nBeta      =  50   % number of temperature steps
        nAnneal    =  20   % number of loops per annealing temp
        
        paramsMap     % parameter name to struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10
        paramsIndices % parameter name to unique integer index
    end

	properties (Dependent)
        fixed
        min
        mean        
        max
        std
        fixedValue
        length
    end
    
    methods %% GET/SET
        function p = get.fixed(this)
            p = this.fixed_;
        end
        function p = get.min(this)
            p = this.min_;
        end
        function p = get.mean(this)
            p = this.mean_;
        end
        function p = get.max(this)
            p = this.max_;
        end
        function p = get.std(this)
            p = this.std_;
        end
        function p = get.fixedValue(this)
            p = this.fixedValue_;
        end
        function p = get.length(this)
            p = this.paramsMap.Count;
        end
    end
    
    methods 
        function this = BayesianParameters(pmap)
            %% BAYESIANPARAMETERS 
            %  Usage:  pmap = containers.Map
            %          pmap('a parameter name') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = BayesianParameters(pmap)
            
            assert(isa(pmap, 'containers.Map') && ~isempty(pmap));
            this.paramsMap = pmap; 
            
            this = this.buildParamsIndices;
            this = this.buildProtectedProperties;
            this.checkParams;
        end
    end
    
    %% PROTECTED
    
    properties (Access = 'protected')
        fixed_
        min_
        mean_
        max_
        std_
        fixedValue_
    end
    
    methods (Access = 'protected')  
        function this = buildParamsIndices(this)
            keys = this.paramsMap.keys;            
            this.paramsIndices = containers.Map;
            for k = 1:length(keys)
                this.paramsIndices(keys{k}) = k;
            end
        end
        function this = buildProtectedProperties(this)
            keys = this.paramsMap.keys;
            this.fixed_      = zeros(length(keys), 1);
            this.min_        = zeros(length(keys), 1);
            this.mean_       = zeros(length(keys), 1);
            this.max_        = zeros(length(keys), 1);
            this.std_        = zeros(length(keys), 1);
            this.fixedValue_ = zeros(length(keys), 1);
            for k = 1:length(keys)
                this.fixed_(k)      = this.paramsMap(keys{k}).fixed;
                this.min_(k)        = this.paramsMap(keys{k}).min;
                this.mean_(k)       = this.paramsMap(keys{k}).mean;
                this.max_(k)        = this.paramsMap(keys{k}).max;
                this.std_(k)        = this.paramsStd( keys{k});
                this.fixedValue_(k) = this.paramsFixedValue(keys{k});
            end
        end
        function s    = paramsStd(this, key)
            s = 0.25*(this.paramsMap(key).max - this.paramsMap(key).min);
        end
        function v    = paramsFixedValue(this, key)
            v = this.paramsMap(key).mean;
        end
        function        checkParams(this)
            min  = this.min_;
            mean = this.mean_;
            max  = this.max_;
            if (any(mean < min) || any(max < mean))
                fprintf('\n');
                fprintf('ERROR:  Prior mean is outside [min max].\n');
                this.printParams;
                error('mlbayesian:parameterOutOfRange', 'BayesianParameters.checkParams');
            end
            if (any(max < min))                
                fprintf('\n');
                fprintf('ERROR:  Prior min > max.\n');
                this.printParams;
                error('mlbayesian:parameterOutOfRange', 'BayesianParameters.checkParams');
            end
            if (any(max == min))
                fprintf('\n');
                fprintf('ERROR:  Parameters will create fault attempting to access mlbayesian.MCMC.lpPopulations(0); \n');
                fprintf('index must be a positive integer or logical.\n');
                fprintf('Fault in mlbayesian.MCMC.replacePoorMembers (line 473).\n');
                fprintf('Try setting prior min != max.\n');
                this.printParams;
                error('mlbayesian:parameterOutOfRange', 'BayesianParameters.checkParams');
            end
        end
        function        printParams(this)
            keys = this.paramsMap.keys;  
            fprintf('\n');
            fprintf('Param. Key \t\tMin \tMean \tMax\n');
            for p = 1:length(keys)
                aMap = this.paramsMap(keys{p});
                fprintf('%s \t\t%g\t%g\t%g\n', keys{p}, aMap.min, aMap.mean, aMap.max);
            end
            fprintf('\n');
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

