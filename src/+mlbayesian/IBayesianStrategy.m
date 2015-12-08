classdef (Abstract) IBayesianStrategy 
	%% IBAYESIANSTRATEGY is the context for AbstractBayesianStrategy 

	%  $Revision$
 	%  was created 23-Nov-2015 17:08:15
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	
    properties (Abstract)
        xLabel
        yLabel
        baseTitle
        detailedTitle
        
        independentData % cells, e.g., times
        dependentData   % cells, e.g., densities = f(time)
        theParameters
        theSolver
        
        expectedBestFitParams
        bestFitParams
        meanParams
        stdParams
        stdOfError
        verbosity
    end
    
    methods (Abstract)
        estimateParameters(this)
        estimateData(this)
        estimateDataFast(this)        
        plot(this)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

