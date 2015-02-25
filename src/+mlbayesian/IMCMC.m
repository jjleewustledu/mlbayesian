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
        NPROPOSALS % number of loops in parameter prob phase
        NPOP       % number of population
        NPOPREP    % number of population to replace
        NBETA      % number of temperature steps
        NANNEAL    % number of loops per annealing temp
        NMOD       % number of models (unused)
    end
       
    properties (Abstract)
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
        
        paramsHist   
        logProbQC   
        stdOfError
        
        nParams
        nSamples
        nProposalsQC
    end

	methods (Abstract)
        runMcmc(this)
        logProbability(this, paramsVec, beta, logProbabilityFlag)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

