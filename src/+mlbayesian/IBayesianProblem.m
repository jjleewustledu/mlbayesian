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
        independentData % numeric
        dependentData   % numeric     
        paramsManager   % mlbayesian.IBayesianParameters
        mcmc            % mlbayesian.IMCMC
        bestFitParams   % numeric
 	end 

	methods (Abstract)
        estimateParameters(this)
  		sumSquaredErrors(this)
        estimateData(this)
        estimateDataFast(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

