classdef IMcmcStrategy 
	%% IMCMCSTRATEGY  

	%  $Revision$
 	%  was created 26-Nov-2015 22:26:04
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	
	properties (Abstract)
        showAnnealing
        showBeta
        showPlots % boolean
        
        dt
        length           % of dependent_data = f(time_interpolants), which must have the same array sizes
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
        annealingSdpar
        annealingInitz
 	end 

	methods (Abstract) 
        finalParams(this)
        finalMeans(this)
        finalStds(this)
        ensureKeyOrdering(this)
        
        runMcmc(this)  
        sumSquaredErrors(this, p)
        adjustParams(this)
        Q(this)
        normalizedQ(this)
        printQNQ(this)
        printBestFit(this)
        printFinalStats(this)
        histParametersDistributions(this)
        histStdOfError(this)
        plotAnnealing(this)
        plotLogProbabilityQC(this)        
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

