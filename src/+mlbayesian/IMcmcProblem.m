classdef (Abstract) IMcmcProblem  
	%% IMCMCPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b).  Copyright 2014 John Joowon Lee. 
 	%% $Id$ 
 	 

	properties (Abstract) 
        baseTitle
        xLabel
        yLabel        
        
        dt
        times            % synonym of independentData
        timeInterpolants % synonym of independentData        
        timeInitial      % independentData(1)
        timeFinal        % independentData(end)
        
        NPROPOSALS % number of loops in parameter prob phase
        NPOP       % number of population
        NPOPREP    % number of population to replace
        NBETA      % number of temperature steps
        NANNEAL    % number of loops per annealing temp        
        
        annealingAvpar
        annealingSdpar
        annealingInitz
        
        showAnnealing
        showBeta
        showPlots % boolean
 	end 

	methods (Abstract) 
        runMcmc(this)
        length(this) % of dependent_data = f(time_interpolants), which must have the same array sizes
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

