classdef IBayesianProblem  
	%% IBAYESIANPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
 	 

	properties (Abstract)         
        independentData % numeric, e.g., times
        dependentData   % numeric, e.g., densities = f(time)
        paramsManager   % mlbayesian.IBayesianParameters
        mcmc            % mlbayesian.IMCMC
        bestFitParams   % numeric
 	end 

	methods (Abstract)
        estimateParameters(this)
  		sumSquaredErrors(this) % merit function
        estimateData(this)     % objective interface for human readability
        estimateDataFast(this) % speed-optimal interface
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

