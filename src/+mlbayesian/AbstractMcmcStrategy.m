classdef AbstractMcmcStrategy < mlbayesian.AbstractBayesianStrategy & mlbayesian.IMcmcStrategy
	%% ABSTRACTMCMCSTRATEGY provides data and behavior for MCMC.
    %  Uses a form of Jeffreys' prior in sumSquaredErrors.

	%  $Revision$
 	%  was created 23-Nov-2015 17:37:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    
    properties    
        cost           = 'sumSquaredResiduals'
        showAnnealing  = true
        showBeta       = true
        showBestFit    = true
        showFinalStats = true
        showPlots      = true
    end
    
    properties (Dependent)        
        annealingAvpar
        annealingInitz
        annealingSdpar        
        bestFitParams
        expectedBestFitParams
        jeffreysPrior
        meanParams
        stdParams
        stdOfError        
              
        nParams
        nProposals % number of proposals for importance sampling, default 100
        nPop       % number of population for annealing/burn-in and proposal/sampling, default 50
        nPopRep    % number of population to replace, default nPop/10
        nBeta      % number of temperature steps, default 50; incr. for more precisions
        nAnneal    % number of loops per annealing temp, default 20; incr. for more precisions
        nSamples   % numel of this.independentData
        nProposalsQC
    end
    
    methods 
        
        %% GET
        
        function a  = get.annealingAvpar(this)
            a = this.theSolver.annealingAvpar; % query what was actually used for MCMC
        end
        function a  = get.annealingInitz(this)
            a = this.theSolver.annealingInitz; % query what was actually used for MCMC
        end
        function a  = get.annealingSdpar(this)
            a = this.theSolver.annealingSdpar; % query what was actually used for MCMC
        end        
        function g  = get.bestFitParams(this)
            assert(~isempty(this.theSolver));
            g = this.theSolver.bestFitParams;
        end
        function g  = get.expectedBestFitParams(this)
            g = cell2mat(this.keysArgs_);
        end
        function j  = get.jeffreysPrior(this)
            j = this.jeffreysPrior_;
        end
        function p  = get.meanParams(this)
            p = this.theSolver.meanParams;
        end
        function g  = get.stdParams(this)
            g = this.theSolver.stdParams;
        end
        function g  = get.stdOfError(this)
            g = this.theSolver.stdOfError;
        end
        
        function n  = get.nParams(this)
            n = this.theParameters.length;
        end
        function n  = get.nProposals(this)
            n = this.theSolver.nProposals; % query what was actually used for MCMC
        end
        function n  = get.nPop(this)
            n = this.theSolver.nPop; % query what was actually used for MCMC
        end
        function n  = get.nPopRep(this)
            n = this.theSolver.nPopRep; % query what was actually used for MCMC
        end
        function n  = get.nBeta(this)
            n = this.theSolver.nBeta; % query what was actually used for MCMC
        end
        function n  = get.nAnneal(this)
            n = this.theSolver.nAnneal; % query what was actually used for MCMC
        end
        function n  = get.nSamples(this)
            n = numel(cell2mat(this.independentData));
        end
        function n  = get.nProposalsQC(this)
            n = this.theSolver.nProposalsQC; % query what was actually used for MCMC
        end
        
        %%
        
        function this = adjustN(this, kind, n)
            assert(ischar(kind));
            assert(isnumeric(n));
            tp = this.theParameters;
            if (isempty(tp))
                tp = mlbayesian.McmcParameters(this.mapParams, numel(cell2mat(this.independentData)));
            end
            tp.(kind) = tp.(kind)*n;
            this.theParameters = tp;
        end
        function ps   = adjustParams(~, ps)
            %% ADJUSTPARAMS:  override as needed for parameter constraints
        end
        function        disp(this)
            builtin('disp', this);
            fprintf('  with mapParams:\n\n');
            this.dispMapParams;
        end
        function        dispMapParams(this)
            keys = this.mapParams.keys;
            vals = this.mapParams.values;
            lenKeys = max(cellfun(@(x) length(x), keys));
            s = '';
            for p = 1:this.mapParams.Count
                s = [s sprintf('    %*s: %s\n', lenKeys, keys{p}, struct2str(vals{p}, 'Punctuation', false))];
            end
            disp(s);
        end
        function        ensureKeyOrdering(this, currentKeys)
            expectedKeys = this.theParameters.keysParams;
            for k = 1:length(expectedKeys)
                assert(strcmp(expectedKeys{k}, currentKeys{k}), ...
                       sprintf('AbstractMcmcStrategy.ensureKeyOrdering:  expected %s but received %s', ...
                       expectedKeys{k}, currentKeys{k}));
            end
        end
        function ed   = estimateData(this)
            keys        = this.theParameters.keysParams;
            finalParams = cellfun(@(x) this.finalParams(x), keys);
            ed          = this.estimateDataFast(finalParams{:});
        end
        function this = estimateParameters(this, varargin)
            ip = inputParser;
            addOptional(ip, 'mapParams', this.mapParams, @(x) isa(x, 'containers.Map'));
            parse(ip, varargin{:});
            
            this = this.runMcmc(ip.Results.mapParams, 'keysToVerify', this.keysParams_);
        end
        function x    = finalMeans(this, key)
            x = this.meanParams(this.theParameters.paramsIndices(key));
        end   
        function x    = finalParams(this, key)
            x = this.bestFitParams(this.theParameters.paramsIndices(key));
        end   
        function x    = finalStds(this, key)
            x = this.stdParams(this.theParameters.paramsIndices(key));
        end      
        function        histParametersDistributions(this) 
            this.theSolver.histParametersDistributions;
        end
        function        histStdOfError(this) 
            this.theSolver.histStdOfError;
        end
        function mdl  = itsModel(this)
            mdl = this.estimateDataFast(this.keysArgs_);
        end
        function lg   = logging(this)
            lg = mlpipeline.Logger(this.fqfileprefix);
            if (isempty(this.summary))
                return
            end
            s = this.summary;
            lg.add('\n%s is working in %s\n', mfilename, pwd);
            lg.add('\begin{alignat*}');
            if (~isempty(this.theSolver))
                lg.add('bestFitParams &-> %s \n', mat2str(s.bestFitParams));
                lg.add('meanParams    &-> %s \n', mat2str(s.meanParams));
                lg.add('stdParams     &-> %s \n', mat2str(s.stdParams));
                lg.add('anneal sdpar  &-> %s \n', mat2str(s.annealingSdpar));
            end
            lg.add('\end{alignat*}');
            lg.add('\n');
        end
        function len  = length(this)
            %% LENGTH
            %  @returns length of dependent_data and of independent_data, both of which must have the same array sizes
            
            for iidx = 1:length(this.independentData)
                len(iidx) = length(this.independentData); %#ok<AGROW>
            end
        end
        function this = makeQuiet(this)
            setenv('VERBOSITY', '0');
            this.showAnnealing  = false;
            this.showBeta       = false;
            this.showBestFit    = false;
            this.showFinalStats = false;
            this.showPlots      = false;
        end
        function nq   = normalizedQ(this)
            aucs = 0;
            for iidx = 1:length(this.dependentData)
                aucs = aucs + sum(abs(this.dependentData{iidx}).^2);
            end
            nq = this.Q/aucs;
        end
        function        plot(~, varargin)
        end
        function        plotAll(this)
            if (~this.showAnnealing)
                this.plotAnnealing;
            end
            if (~this.showPlots)
                this.histParametersDistributions;
                this.plotLogProbabilityQC;
                this.histStdOfError;
            end
            this.plot;
        end
        function        plotAnnealing(this) 
            this.theSolver.plotAnnealing;
        end
        function        plotLogProbabilityQC(this) 
            this.theSolver.plotLogProbabilityQC;
        end
        function        plotParVars(~, varargin)
        end
        function        printBestFit(this)
            this.theSolver.printBestFit;
        end
        function        printFinalStats(this)
            this.theSolver.printFinalStats;
        end
        function        printQNQ(this)
            fprintf('FINAL STATS Q               %g\n', this.Q);
            fprintf('FINAL STATS Q normalized    %g\n', this.normalizedQ);
        end
        function q    = Q(this)
            q = this.sumSquaredErrors(this.bestFitParams);
        end        
        function this = runMcmc(this, varargin)
            %% RUNMCMC should be run from within method estimateParameters, implemented as described below.
            %  @params mapParams is a containers.Map with parameter keys & values
            %  @params named keysParams is a cell array enumerating param keys to be checked for equivalence with
            %  mapParams.keys
            %  @params named mcmcParameters may be assigned in mlbayesian.AbstractBayesianStrategy.ctor,
            %  by direct assignment of this.theParameters or here; reasonable defaults are used otherwise.
            %  Usage:  map = containers.Map
            %          map('some_param') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = this.runMcmc(map)
            
            import mlbayesian.*;  
            ip = inputParser;
            addOptional( ip, 'mapParams', this.mapParams, @(x) isa(x, 'containers.Map'));
            addParameter(ip, 'keysParams', {}, @iscell);
            addParameter(ip, 'mcmcParameters', this.theParameters_, @(x) isa(x, 'mlbayesian.IMcmcParameters'));
            parse(ip, varargin{:});                    
            if (~isempty(ip.Results.keysParams))
                this.ensureKeyOrdering(ip.Results.keysParams);
            end            
            this.theParameters_ = ip.Results.mcmcParameters;  
            if (isempty(this.theParameters_))
                this.theParameters_ = McmcParameters( ...
                    ip.Results.mapParams, numel(cell2mat(this.independentData)), ...
                    'pkeys', ip.Results.keysParams);
            end
            
            this.theSolver       = McmcCellular(this);
            [~,~,this.theSolver] = this.theSolver.runMcmc;     
            keys_ = ip.Results.keysParams;
            for k = 1:length(keys_) 
                this.(keys_{k}) = this.finalParams(keys_{k});
            end
            this.printQNQ;
        end
        function this = simulateItsMcmc(this)
        end
        function sse  = sumSquaredErrors(this, p)
            %% SUMSQUAREDERRORS 
            
            switch (this.cost)
                case 'sumSquaredResiduals'
                    sse = this.sumSquaredResiduals(p);
                case 'sumSquaredResidualsJeffreys'
                    sse = this.sumSquaredResidualsJeffreys(p);
                otherwise
                    error('mlbayesian:unsupportedSwitchCase', 'AbstractMcmcStrategy.sumSquaredErrors');
            end
        end
        function this = updateSummary(this)
            s.class = class(this);
            s.datestr = datestr(now, 30);
            if (~isempty(this.theSolver))
                s.bestFitParams  = this.bestFitParams;
                s.meanParams     = this.meanParams;
                s.stdParams      = this.stdParams;
                s.annealingSdpar = this.annealingSdpar;
            end
            this.summary = s;
        end
        
 		function this = AbstractMcmcStrategy(varargin) 
 			%% ABSTRACTMCMCSTRATEGY 
 			%  Usage:  this = AbstractMcmcStrategy() 

 			this = this@mlbayesian.AbstractBayesianStrategy(varargin{:});
        end 
        
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        jeffreysPrior_
        keysArgs_
        keysParams_
        mapParams_
    end
    
    methods (Access = protected)
        function this = buildJeffreysPrior(this)
            %% JEFFREYSPRIOR
            %  Cf. Gregory, Bayesian Logical Data Analysis for the Physical Sciences, sec. 3.7.1.
            
            p = cell(this.independentData);
            for iidx = 1:length(p)
                t = this.independentData{iidx};
                for it = 1:length(t)
                    if (abs(t(it)) < eps)
                        t(it) = min(t(t > eps));
                    end
                end
                p{iidx} = 1./(t*log(t(end)/t(1)));
                p{iidx} = p{iidx}/sum(p{iidx});
            end
            this.jeffreysPrior_ = p;
            this.cost = 'sumSquaredResidualsJeffreys';
        end
        function this = buildJeffreysPriorNonuniform(this)
            %% JEFFREYSPRIOR
            %  Cf. Gregory, Bayesian Logical Data Analysis for the Physical Sciences, sec. 3.7.1.
            
            p = cell(this.independentData);
            for iidx = 1:length(p)
                t = this.independentData{iidx};
                for it = 1:length(t)
                    if (abs(t(it)) < eps)
                        t(it) = min(t(t > eps));
                    end
                end
                taus_ = t(2:end) - t(1:end-1);
                taus_ = [taus_ taus_(end)]; %#ok<AGROW>
                p{iidx} = 1./(t*log(t(end)/t(1)));
                p{iidx} = p{iidx}.*taus_/sum(p{iidx}.*taus_);
            end
            this.jeffreysPrior_ = p;
            this.cost = 'sumSquaredResidualsJeffreys';
        end
        function s    = keysParamsForSprintf(this, argsv)
            assert(length(argsv) == length(this.keyParams_));
            s = '';
            for a = 1:length(argsv)-1
                s = [s this.keysParams_{a} ' %g, '];
            end
            s = [s this.keysParams_{end} ' %g'];
        end
        function        plotParArgs(this, par, args, vars)
            assert(lstrfind(par, properties(this)));
            assert(iscell(args));
            assert(isnumeric(vars));
            figure
            hold on
            mdl = this.itsModel;
            for v = 1:length(args)
                argsv = args{v};
                plot(0:length(mdl)-1, ...
                     mlbayesian.AbstractMcmcStrategy.model(argsv{:}, this.times{1}));
            end
            plot(0:length(mdl)-1, this.dependentData{1}, 'o', 'LineWidth', 2);
            title(sprintf(this.keysParamsForSprintf(argsv), argsv{:}));
            legend(['bayes' ...
                    cellfun(@(x) sprintf('%s = %g', par, x), num2cell(vars), 'UniformOutput', false)]);
            xlabel('time sampling index');
            ylabel(this.yLabel);
        end        
        function sse  = sumSquaredResiduals(this, p)
            %% SUMSQUAREDRESIDUALS returns the sum-of-square residuals for all cells of this.dependentData and 
            %  corresponding this.estimateDataFast.  
            
            p   = num2cell(p);
            sse = 0;
            edf = this.estimateDataFast(p{:});
            for iidx = 1:length(this.dependentData)
                sse = sse + ...
                      sum( (this.dependentData{iidx} - edf{iidx}).^2 );
            end
            if (sse < 10*eps)
                sse = sse + (1 + rand(1))*10*eps; 
            end
        end        
        function sse  = sumSquaredResidualsJeffreys(this, p)
            %% SUMSQUAREDRESIDUALSJEFFREYS returns the sum-of-square residuals for all cells of this.dependentData and 
            %  corresponding this.estimateDataFast.  This 
            %  implementation weights the sum-of-square residuals with Jeffrey's prior according to this.independentData.
            
            p   = num2cell(p);
            sse = 0;
            edf = this.estimateDataFast(p{:});
            for iidx = 1:length(this.dependentData)
                sse = sse + ...
                      sum( (this.dependentData{iidx} - edf{iidx}).^2 .* ...
                            this.jeffreysPrior_{iidx} );
            end
            if (sse < 10*eps)
                sse = sse + (1 + rand(1))*10*eps; 
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

