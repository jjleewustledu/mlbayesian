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
        showAnnealing = false
        showBeta      = false
        showPlots     = false
    end
    
    properties (Dependent)
        dt
        length           % of dependent_data and of indepdendent_data, both of which must have the same array sizes
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
    
    methods %% GET
        function t  = get.dt(this)            
            for iidx = 1:length(this.independentData)
                t{iidx} = min(this.taus{iidx});
            end
        end
        function le = get.length(this)
            for iidx = 1:length(this.independentData)
                le(iidx) = length(this.independentData);
            end
        end
        function t  = get.taus(this)            
            for iidx = 1:length(this.independentData)
                t{idx} = this.times{idx}(2:end) - this.times{idx}(1:end-1);
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
        function n = get.nSamples(this)
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
    end

	methods
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
            storedKeys = this.theParameters.paramsMap.keys;
            for k = 1:length(storedKeys)
                assert(strcmp(storedKeys{k}, currentKeys{k}), ...
                       sprintf('AbstractMcmcStrategy.ensureKeyOrdering:  expected %s but received %s', ...
                       storedKeys{k}, currentKeys{k}));
            end
        end
        function ed   = estimateData(this)
            keys        = this.theParameters.paramsMap.keys;
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
        function nq   = normalizedQ(this)
            aucs = 0;
            for iidx = 1:length(this.dependentData)
                aucs = aucs + sum(abs(this.dependentData{iidx}).^2);
            end
            nq = this.Q/aucs;
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
            %  @params paramsMap is a containers.Map with parameter keys & values
            %  @params paramsKeys is a cell array enumerating param keys to be checked for equivalence with
            %  paramsMap.keys
            %  Usage:  map = containers.Map
            %          map('some_param') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = this.runMcmc(map)
            
            ip = inputParser;
            addRequired( ip, 'paramsMap',        @(x) isa(x, 'containers.Map'));
            addParameter(ip, 'keysToVerify', {}, @iscell);
            parse(ip, varargin{:});
            
            import mlbayesian.*;
            this.theParameters   = McmcParameters(ip.Results.paramsMap, numel(cell2mat(this.independentData)));    
            if (~isempty(ip.Results.keysToVerify))
                this.ensureKeyOrdering(ip.Results.keysToVerify);
            end
            this.theSolver       = McmcCellular(this);
            [~,~,this.theSolver] = this.theSolver.runMcmc;     
            this.printQNQ;
            keys = ip.Results.paramsMap.keys;
            for k = 1:length(keys) %#ok<CPROPLC>
                this.(keys{k}) = this.finalParams(keys{k});
            end
        end
        function sse  = sumSquaredErrors(this, p)
            %% SUMSQUAREDERRORS returns the sum-of-squared-errors summed over the cells of this.dependentData and 
            %  corresponding this.estiamteDataFast. 
            
            p   = num2cell(p);
            sse = 0;
            edf = this.estimateDataFast(p{:});
            for iidx = 1:length(this.dependentData) %#ok<CPROPLC>
                sse = sse + ...
                      sum(abs(this.dependentData{iidx} - edf{iidx}).^2) ./ sum(abs(this.dependentData{iidx}).^2);
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

