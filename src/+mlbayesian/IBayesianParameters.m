classdef IBayesianParameters  
	%% IBAYESIANPARAMETERS is used by implementations of mlbayesian.IMCMC  

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a).  Copyright 2014 John Joowon Lee. 
 	%  $Id$  	 

	properties (Abstract)
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

