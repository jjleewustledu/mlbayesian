classdef (Abstract) IMcmcParameters 
	%% IMCMCPARAMETERS  

	%  $Revision$
 	%  was created 23-Nov-2015 18:22:18
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	
	properties (Abstract)
        nProposals % number of loops in parameter prob phase
        nPop       % number of population
        nBeta      % number of temperature steps
        nAnneal    % number of loops per annealing temp
        nSamples
        
        paramsMap     % parameter name to struct that defines values for fixed, min, mean, max, std, fixedValue
        paramsIndices % parameter name to unique integer index
        
        fixed
        min
        mean
        max
        std
        fixedValue
        length
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

