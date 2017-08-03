classdef PolynomialMcmcCellular < mlbayesian.AbstractMcmcStrategy
	%% POLYNOMIALMCMCCELLULAR is for testing McmcCellular implementations; McmcCellular requires access to an AbstractMcmcStrategy

	%  $Revision$
 	%  was created 23-Nov-2015 17:37:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
    
    
    properties
        c0 = 1
        c1 = 1
        c2 = 1
        
        xLabel = 't'
        yLabel = '\text{poly}^3(t)'
    end
    
    properties (Dependent)
        detailedTitle
        mapParams        
    end
    
    methods (Static)        
        function mdl  = model(c0, c1, c2, t)
            mdl = c0 + c1*t + c2*t.^2;
        end
    end
    
	methods
        
        %% GET
        
        function dt = get.detailedTitle(this)
            dt = sprintf('%s:\nc0 %g, c1 %g, c2 %g', ...
                         this.baseTitle, ...
                         this.c0, this.c1, this.c2);
        end
        function m  = get.mapParams(this)
            m = containers.Map;
            m('c0') = struct('fixed', 0, 'min', eps, 'mean', this.c0, 'max',  10*this.c0);
            m('c1') = struct('fixed', 0, 'min', eps, 'mean', this.c1, 'max',  10*this.c1);
            m('c2') = struct('fixed', 0, 'min', eps, 'mean', this.c2, 'max',  10*this.c2);
        end
        
        %%
        
        function this = PolynomialMcmcCellular(varargin)
            this = this@mlbayesian.AbstractMcmcStrategy(varargin{:});
        end
        
        function this = estimateParameters(this)
            this = this.runMcmc(this.mapParams);
        end   
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.bestFitParams(1), this.bestFitParams(2), this.bestFitParams(3));
        end
        function ed   = estimateDataFast(this, c0, c1, c2)
            ed{1} = this.model(c0, c1, c2, this.independentData{1});
        end  
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

