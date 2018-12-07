classdef (Abstract) ISolver 
	%% ISOLVER  

	%  $Revision$
 	%  was created 12-Dec-2017 17:36:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlanalysis/src/+mlanalysis.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties (Abstract)		
 		useSynthetic % used stored parameters to generate synthetic data
        isfinished
        model
 	end

	methods (Abstract)
        diagnose(this)
        this = estimateParameters(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

