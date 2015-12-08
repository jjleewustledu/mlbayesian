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
        paramsBetas
        paramsPopulations  
        paramsSigmas
        paramsHist 
        annealingAvpar
        annealingSdpar
        annealingInitz
        bestFitParams
        meanParams
        stdParams
        stdOfError
        
        lpQC
        lpBetas
        lpPopulations
        lpFinal
        
        parameters        
        nParams
        nProposals % number of loops in parameter prob phase
        nPop       % number of population        
        nPopRep    % number of population to replace
        nBeta      % number of temperature steps
        nAnneal    % number of loops per annealing temp
        nSamples
        nProposalsQC
        showAnnealing
        showBeta
        showPlots
        verbosity
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

