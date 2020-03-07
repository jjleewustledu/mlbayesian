classdef IBretthorstModel < mlbayesian.IModel
	%% IBRETTHORSTMODEL  

	%  $Revision$
 	%  was created 15-Sep-2018 00:42:49 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties (Abstract)
 		independentData
        dependentData        
        useSynthetic % let initial parameters form synthetic data 		
 	end

	methods (Abstract)
        this = estimateData(this)
        Q    = objectiveFunc(this)        
        ps   = modelParameters(this)    
        sps  = modelStdParameters(this)
        ps   = solverParameters(this)
        this = updateModel(this, solvr)
        
        plot(this)
        writetable(this)		  
  	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

