classdef SumSquaredWeightedResiduals < handle & mlbayesian.SumSquaredResiduals
	%% SUMSQUAREDWEIGHTEDRESIDUALS  

	%  $Revision$
 	%  was created 29-Aug-2018 00:02:33 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

	methods        
        function e = sumSquaredErrors(this, p)
            %% SUMSQUAREDERRORS returns the sum-of-square residuals weighted by 1/sigmas of a Poisson distribution for 
            %  noise for all cells of this.client_.dependentData and this.client_.estimateDataFast. 
            p = num2cell(p);
            e = 0;
            edf = this.client_.estimateDataFast(p{:});
            for idd = 1:length(this.client_.dependentData)
                summand = abs(this.client_.dependentData{idd} - edf{idd}).^2 ./ ...
                              this.client_.dependentData{idd};
                e = e + sum(summand(isfinite(summand)));
            end
            if (e < 10*eps)
                e = e + (1 + rand(1))*10*eps; 
            end
        end		
        
 		function this = SumSquaredWeightedResiduals(varargin)
 			%% SUMSQUAREDWEIGHTEDRESIDUALS
 			%  @param .

 			this = this@mlbayesian.SumSquaredResiduals(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

