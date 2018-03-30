classdef IMcmcSolver < mlanalysis.ISolver
	%% IMCMCSOLVER  

	%  $Revision$
 	%  was created 02-Jan-2018 20:19:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

	methods (Abstract)		
        plotAnnealing(this)
        plotLogProbabilityQC(this)
        plotParameterCovariances(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy    
    
 end

