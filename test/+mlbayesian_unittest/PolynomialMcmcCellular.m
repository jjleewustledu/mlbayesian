classdef PolynomialMcmcCellular < mlbayesian.AbstractMcmcStrategy
	%% POLYNOMIALMCMCCELLULAR is for testing McmcCellular implementations; McmcCellular requires access to an AbstractMcmcStrategy

	%  $Revision$
 	%  was created 23-Nov-2015 17:37:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
    
    methods (Static)        
        function mdl  = model(varargin)
            mdl = [];
        end
    end
    
	methods
        function this = PolynomialMcmcCellular(indDat, depDat)
            this = this@mlbayesian.AbstractMcmcStrategy(indDat, depDat);
        end
        
        function this = estimateParameters(this)  
            mapParams = containers.Map;
            mapParams('c0') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            mapParams('c1') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            mapParams('c2') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            this = this.runMcmc(mapParams);
        end   
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.bestFitParams(1), this.bestFitParams(2), this.bestFitParams(3));
        end
        function ed   = estimateDataFast(this, c0, c1, c2)
            ed{1} = c0 + c1*this.independentData{1} + c2*this.independentData{1}.^2;
        end  
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

