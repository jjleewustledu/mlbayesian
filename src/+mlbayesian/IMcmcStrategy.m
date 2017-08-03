classdef IMcmcStrategy 
	%% IMCMCSTRATEGY  

	%  $Revision$
 	%  was created 26-Nov-2015 22:26:04
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  Copyright 2015 John Joowon Lee. 
 	
	properties (Abstract)        
        cost
        mapParams
        showAnnealing
        showBeta
        showBestFit
        showFinalStats
        showPlots        
        
        annealingAvpar
        annealingInitz
        annealingSdpar        
        bestFitParams
        expectedBestFitParams
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

	methods (Abstract) 
        this = adjustN(this, kind, n)
        ps   = adjustParams(~, ps)
               disp(this)
               dispMapParams(this)
               ensureKeyOrdering(this, currentKeys)
        ed   = estimateData(this)
        edf  = estimateDataFast(this)
        this = estimateParameters(this, varargin)
        
        x    = finalMeans(this)
        x    = finalParams(this)
        x    = finalStds(this)               
               histParametersDistributions(this)
               histStdOfError(this)
        mdl  = itsModel(this)
        len  = length(this)
        nq   = normalizedQ(this)
               plot(~, varargin)
               plotAll(this)
               plotAnnealing(this)
               plotLogProbabilityQC(this)
               plotParVars(~, varargin)
               printBestFit(this)
               printFinalStats(this)
               printQNQ(this)
        q    = Q(this)  
        this = runMcmc(this)  
        this = simulateItsMcmc(this)
        sse  = sumSquaredErrors(this, p)     
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

