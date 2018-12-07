classdef McmcParameters < mlbayesian.IMcmcParameters 
	%% MCMCPARAMETERS  

	%  $Revision$
 	%  was created 23-Nov-2015 18:22:26
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	
    properties % must be in heap memory for speed
        nProposals = 100 % number of proposals for importance sampling, ~100
        nPop       =  50 % number of population for annealing/burn-in and proposal/sampling, ~50
        nBeta      =  50 % number of temperature steps, ~50
        nAnneal    =  20 % number of loops per annealing temp, ~20
        nSamples         % numel of independentData     
        
        fixed
        fixedValue
        min_
        mean_        
        max_
        std_
    end

	properties (Dependent)        
        indicesParams % parameter name to unique integer index
        keysParams    % parameter keys with cell ordering
        mapParams     % parameter name to mapped struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
        paramsIndices % synonym for legacy support
        paramsKeys    % "
        paramsMap     % "
    end
    
    methods 
        
        %% GET
                
        function g    = get.indicesParams(this)
            g = this.indicesParams_;
        end
        function g    = get.keysParams(this)
            g = this.keysParams_;
        end
        function g    = get.mapParams(this)
            g = this.mapParams_;
        end
        function g    = get.paramsIndices(this)
            g = this.indicesParams_;
        end
        function g    = get.paramsKeys(this)
            g = this.keysParams_;
        end
        function g    = get.paramsMap(this)
            g = this.mapParams_;
        end

        %%
        
        function len  = length(this)
            len = this.mapParams_.Count;
        end
        
        function this = McmcParameters(varargin)
            %% MCMCPARAMETERS 
            %  Usage:  pmap = containers.Map
            %          pmap('a parameter name') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = McmcParameters(pmap)
            
            ip = inputParser;
            addRequired( ip, 'pmap',     @(x) isa(x, 'containers.Map') && ~isempty(x));
            addRequired( ip, 'nsamples', @isnumeric);
            addParameter(ip, 'pkeys', {}, @iscell);
            parse(ip, varargin{:});
            
            this.mapParams_ = ip.Results.pmap; 
            this.nSamples = ip.Results.nsamples;
            this.keysParams_ = this.mapParams_.keys; % ordering of keys will be alphanumeric, caps first
            if (~isempty(ip.Results.pkeys))
                this.keysParams_ = ip.Results.pkeys;
            end
            
            this = this.buildParamsIndices;
            this = this.buildProtectedProperties;
            this.checkParams;
        end
    end
        
    %% PRIVATE
    
    properties (Access = 'private')
        indicesParams_
        keysParams_
        mapParams_
    end
    
    methods (Access = 'private')  
        function this = buildParamsIndices(this)
            keys = this.keysParams_;
            this.indicesParams_ = containers.Map;
            for k = 1:length(keys)
                this.indicesParams_(keys{k}) = k;
            end
        end
        function this = buildProtectedProperties(this)
            keys = this.keysParams_;
            this.fixed      = false(length(keys), 1);
            this.min_        = zeros(length(keys), 1);
            this.mean_       = zeros(length(keys), 1);
            this.max_        = zeros(length(keys), 1);
            this.std_        = zeros(length(keys), 1);
            this.fixedValue = zeros(length(keys), 1);
            for k = 1:length(keys)
                this.fixed(k)      = logical(this.mapParams_(keys{k}).fixed);
                this.min_(k)        = this.mapParams_(keys{k}).min;
                this.mean_(k)       = this.mapParams_(keys{k}).mean;
                this.max_(k)        = this.mapParams_(keys{k}).max;
                this.std_(k)        = this.paramsStd( keys{k});
                this.fixedValue(k) = this.paramsFixedValue(keys{k});
            end
        end
        function s    = paramsStd(this, key)
            s = 0.25*(this.mapParams_(key).max - this.mapParams_(key).min);
        end
        function v    = paramsFixedValue(this, key)
            v = this.mapParams_(key).mean;
        end
        function        checkParams(this)
            min__  = this.min_;
            mean__ = this.mean_;
            max__  = this.max_;
            if (any(mean__ < min__) || any(max__ < mean__))
                fprintf('\n');
                fprintf('ERROR:  Prior mean is outside [min max].\n');
                this.printParams;
                error('mlbayesian:parameterOutOfRange', 'McmcParameters.checkParams');
            end
            if (any(max__ < min__))                
                fprintf('\n');
                fprintf('ERROR:  Prior min > max.\n');
                this.printParams;
                error('mlbayesian:parameterOutOfRange', 'McmcParameters.checkParams');
            end
            if (any(max__ == min__))
                fprintf('\n');
                fprintf('ERROR:  Parameters will create fault attempting to access mlbayesian.MCMC.lpPopulations(0); \n');
                fprintf('index must be a positive integer or logical.\n');
                fprintf('Fault in mlbayesian.MCMC.replacePoorMembers (line 473).\n');
                fprintf('Try setting prior min != max.\n');
                this.printParams;
                error('mlbayesian:parameterOutOfRange', 'McmcParameters.checkParams');
            end
        end
        function        printParams(this)
            keys = this.mapParams_.keys;  
            fprintf('\n');
            fprintf('Param. Key \t\tMin \tMean \tMax\n');
            for p = 1:length(keys)
                aMap = this.mapParams_(keys{p});
                fprintf('%s \t\t%g\t%g\t%g\n', keys{p}, aMap.min, aMap.mean, aMap.max);
            end
            fprintf('\n');
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

