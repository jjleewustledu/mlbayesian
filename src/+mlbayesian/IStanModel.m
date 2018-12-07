classdef IStanModel < mlbayesian.IModel
	%% ISTANMODEL  

	%  $Revision$
 	%  was created 15-Sep-2018 00:42:27 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		modelCode
        data
        fit
        iter
        chains
        verbose
 	end

	methods 
        extract(this)
        print(this)
        fit = stan(this)
        traceplot(this)
  	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

