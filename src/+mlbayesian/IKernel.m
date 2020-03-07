classdef (Abstract) IKernel 
	%% IKERNEL  

	%  $Revision$
 	%  was created 14-Sep-2018 16:53:56 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties (Abstract)
        independentData
 		dependentData        
 	end

	methods (Abstract, Static)
		  this = main(varargin)
          ps   = adjustParams(ps)
    end 
    
    methods (Abstract)
        ed = estimateData(this, varargin)
        Q  = objectiveFunc(this, varargin)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

