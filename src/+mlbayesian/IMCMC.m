classdef IMCMC  
	%% IMCMC   
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
	properties (Abstract)
        nProposals % number of loops in parameter prob phase
        nPop       % number of population
        nPopRep    % number of population to replace
        nBeta      % number of temperature steps
        nAnneal    % number of loops per annealing temp

        dependentData
        paramsData
        paramsBetas
        paramsPopulations  
        paramsSigmas
        annealingAvpar
        annealingSdpar
        annealingInitz
        bestFitParams
        
        lpBetas
        lpPopulations
        lpFinal
        
        paramsHist   
        logProbQC   
        stdOfError
        
        nParams
        nSamples
        nProposalsQC
        showAnnealing
        showBeta
        showPlots
    end

	methods (Abstract)
        runMcmc(this)
        logProbability(this, paramsVec, beta, logProbabilityFlag)
        
        printBestFit(this)
        printFinalStats(this)
        histParametersDistributions(this)
        histStdOfError(this)
        plotAnnealing(this)
        plotLogProbabilityQC(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

