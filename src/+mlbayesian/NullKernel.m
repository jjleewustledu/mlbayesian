classdef NullKernel 
	%% NULLKERNEL  

	%  $Revision$
 	%  was created 11-Jan-2018 20:30:20 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        independentData
 		dependentData
 	end

    methods (Static)
        function this = main(varargin)
            this = mlbayesian.NullKernel(varargin{:});
        end
        function ps   = adjustParams(ps)
        end
        function t    = nullfunc(t)
        end
    end 

	methods 
        
        function ed = estimateData(this, varargin)
            %% ESTIMATEDATA
            %  @returns ed, the estimated data isomorphic to this.dependentData.
            
            ed = this.nullfunc(this.independentData);
        end
        function Q  = objectiveFunc(this, varargin)
            %% OBJECTIVEFUNC returns the sum-of-square residuals for all cells of this.dependentData and corresponding
            %  this.estimateDataFast.
            
            Q = sum( (this.dependentData - this.estimateData(varargin{:})).^2 );
            if (Q < 10*eps)
                Q = Q + (1 + rand(1))*10*eps; 
            end
        end
		  
 		function this = NullKernel(idata, ddata)
 			%% NULLKERNEL
            
            assert(isnumeric(idata));
            assert(isnumeric(ddata));

 			this.independentData = idata;	
 			this.dependentData   = ddata;
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

