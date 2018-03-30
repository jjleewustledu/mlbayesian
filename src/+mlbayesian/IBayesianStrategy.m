classdef (Abstract) IBayesianStrategy 
	%% IBAYESIANSTRATEGY is the context for AbstractBayesianStrategy 

	%  $Revision$
 	%  was created 23-Nov-2015 17:08:15
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  Copyright 2015 John Joowon Lee. 
 	
    properties (Abstract)
        baseTitle
        dependentData % cells, e.g., densities = f(time)
        dependentDataInterp % cells, interpolated by this.NyquistFreq
        dt % scalar := min(min(cell2mat(this.independentDeltas)))/this.NyquistFreqFactor
        independentData % cells, e.g., times
        independentDataInterp % cells, interpolated by this.NyquistFreq
        independentDeltas % repeats independentData(end) - independentData(end-1) so that 
                          % length(independentDeltas) = length(independentData).
        NyquistFreqFactor
        sessionData
        theParameters
        theSolver
        verbosity
        
        taus             % synonym of independentDeltas
        times            % synonym of independentData
        timeFinal        % independentData(end)
        timeInitial      % independentData(1)
        timeInterpolants % synonym of independentDataInterp
    end
    
    methods (Static, Abstract)
        load
    end
    
    methods (Abstract)
        save(this)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

