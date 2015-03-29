classdef IBayesianProblem  
	%% IBAYESIANPROBLEM is used by IMCMC, MCMC; most properties and some methods are requested by MCMC   

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
        estimateData(this)     % objective interface for human readability
        estimateDataFast(this) % speed-optimal interface
  		sumSquaredErrors(this) % merit function
        adjustParams(this)
        Q(this)
        normalizedQ(this)
        finalParams(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

