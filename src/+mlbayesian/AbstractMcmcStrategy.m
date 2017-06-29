classdef AbstractMcmcStrategy < mlbayesian.AbstractBayesianStrategy & mlbayesian.IMcmcStrategy
	%% ABSTRACTMCMCSTRATEGY  

	%  $Revision$
 	%  was created 23-Nov-2015 17:37:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties (Abstract)
        mapParams
    end
    
    properties  
        jeffreysPrior
        showAnnealing  = true
        showBeta       = true
        showBestFit    = true
        showFinalStats = true
        showPlots      = true
    end
    
    properties (Dependent)
        dt
        taus             % times(2) - times(1), times(3) - times(2), ...
        times            % synonym of independentData
        timeFinal        % independentData(end)
        timeInitial      % independentData(1)
        timeInterpolants % timeInitial:dt:timeFinal
              
        nParams
        nProposals % number of loops in parameter prob phase  
        nPop       % number of population
        nPopRep    % number of population to replace 
        nBeta      % number of temperature steps
        nAnneal    % number of loops per annealing temp   
        nSamples
        nProposalsQC
        
        annealingAvpar
        annealingInitz
        annealingSdpar
    end
    
    methods 
        
        %% GET
        
        function t  = get.dt(this)
            for iidx = 1:length(this.independentData)
                t{iidx} = min(this.taus{iidx});
            end
        end
        function t  = get.taus(this)
            for iidx = 1:length(this.independentData)
                t{idx} = this.times{idx}(2:end) - this.times{idx}(1:end-1); %#ok<AGROW>
            end
        end
        function t  = get.times(this)
            t = this.independentData;
        end
        function tf = get.timeFinal(this)
            for iidx = 1:length(this.independentData)
                tf(iidx) = this.independentData{iidx}(end);
            end
        end
        function tf = get.timeInitial(this) 
            for iidx = 1:length(this.independentData)
                tf(iidx) = this.independentData{iidx}(1);
            end
        end
        function ti = get.timeInterpolants(this)
            for iidx = 1:length(this.independentData)
                ti{iidx} = this.timeInitial{iidx}:this.dt{iidx}:this.timeFinal{iidx};
            end
        end
        
        function n  = get.nParams(this)
            n = this.theParameters.length;
        end
        function n  = get.nProposals(this)
            n = this.theSolver.nProposals;
        end
        function n  = get.nPop(this)
            n = this.theSolver.nPop;
        end
        function n  = get.nPopRep(this)
            n = this.theSolver.nPopRep;
        end
        function n  = get.nBeta(this)
            n = this.theSolver.nBeta;
        end
        function n  = get.nAnneal(this)
            n = this.theSolver.nAnneal;
        end
        function n  = get.nSamples(this)
            n = numel(cell2mat(this.independentData));
        end
        function n  = get.nProposalsQC(this)
            n = this.theSolver.nProposalsQC;
        end
        
        function a  = get.annealingAvpar(this)
            a = this.theSolver.annealingAvpar;
        end
        function a  = get.annealingInitz(this)
            a = this.theSolver.annealingInitz;
        end
        function a  = get.annealingSdpar(this)
            a = this.theSolver.annealingSdpar;
        end
        
        %%
        
 		function this = AbstractMcmcStrategy(varargin) 
 			%% ABSTRACTMCMCSTRATEGY 
 			%  Usage:  this = AbstractMcmcStrategy() 

 			this = this@mlbayesian.AbstractBayesianStrategy(varargin{:}); 
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
        function x    = finalParams(this, key)
            x = this.bestFitParams(this.theParameters.paramsIndices(key));
        end   
        function x    = finalMeans(this, key)
            x = this.meanParams(this.theParameters.paramsIndices(key));
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
        function len  = length(this)
            %% LENGTH
            %  @returns length of dependent_data and of indepdendent_data, both of which must have the same array sizes
            
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
            %  @params paramsKeys is a cell array enumerating param keys to be checked for equivalence with
            %  mapParams.keys
            %  Usage:  map = containers.Map
            %          map('some_param') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = this.runMcmc(map)
            
            ip = inputParser;
            addRequired( ip, 'mapParams',        @(x) isa(x, 'containers.Map'));
            addParameter(ip, 'keysToVerify', {}, @iscell);
            parse(ip, varargin{:});
            
            import mlbayesian.*;
            this.theParameters   = McmcParameters( ...
                ip.Results.mapParams, numel(cell2mat(this.independentData)), ...
                'pkeys', ip.Results.keysToVerify);    
            if (~isempty(ip.Results.keysToVerify))
                this.ensureKeyOrdering(ip.Results.keysToVerify);
            end
            this.theSolver       = McmcCellular(this);
            [~,~,this.theSolver] = this.theSolver.runMcmc;     
            this.printQNQ;
            keys_ = ip.Results.keysToVerify;
            for k = 1:length(keys_) 
                this.(keys_{k}) = this.finalParams(keys_{k});
            end
        end
        function sse  = sumSquaredErrors(this, p)
            %% SUMSQUAREDERRORS returns the sum-of-square residuals for all cells of this.dependentData and 
            %  corresponding this.estimateDataFast.  Compared to AbstractMcmcStrategy.sumSquaredErrors, this 
            %  overriding implementation weights of the log-likelihood with Jeffrey's prior according to this.independentData.
            %  See also:  mlbayesian.AbstractMcmcStrategy.sumSquaredErrors, 
            %             mlkinetics.AbstractKinetics.jeffreysPrior.
            
            p   = num2cell(p);
            sse = 0;
            edf = this.estimateDataFast(p{:});
            for iidx = 1:length(this.dependentData)
                %sigma2 = this.dependentData{iidx};
                %sigma2(sigma2 < eps) = 1;
                sse = sse + ...
                      sum( (this.dependentData{iidx} - edf{iidx}).^2 .* ...
                            this.jeffreysPrior{iidx} );
                          % ./ sigma2 );
            end
            if (sse < 10*eps)
                sse = sse + (1 + rand(1))*10*eps; 
            end
        end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
    end
    
    methods (Access = protected)
        function p = buildJeffreysPrior(this)
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
        end
        function p = buildJeffreysPrior__(this)
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
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

