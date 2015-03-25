classdef IMcmcProblem  
	%% IMCMCPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 

	properties (Abstract)  
        baseTitle
        xLabel
        yLabel        
        
        length % of dependent_data = f(time_interpolants), which must have the same array sizes
        timeInterpolants
        timeFinal
        
        NPROPOSALS % number of loops in parameter prob phase
        NPOP       % number of population
        NPOPREP    % number of population to replace
        NBETA      % number of temperature steps
        NANNEAL    % number of loops per annealing temp        
        
        annealingAvpar
        annealingSdpar
        annealingInitz
        
        showPlots % boolean
 	end 

	methods (Abstract)
        runMcmc(this)
        logProbability(this, paramsVec, beta, logProbabilityFlag) %  lPFlag: -1  return Q
                                                                  %           0  return log of normalized Q
                                                                  %           1  return with beta and priors
        printBestFit(this)
        printFinalStats(this)
        histParametersDistributions(this)
        histStdOfError(this)
        plotAnnealing(this)
        plotLogProbabilityQC(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

