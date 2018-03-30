classdef PolynomialsKernel < mlbayesian.NullKernel
	%% POLYNOMIALSKERNEL needs small data structures that will all fit in cpu caches. 

	%  $Revision$
 	%  was created 24-Dec-2017 17:52:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.

    methods (Static)
        function this = main(varargin)
            this = mlbayesian.PolynomialsKernel(varargin{:});
        end
        function y    = polynomial(as, t)
            %% POLYNOMIAL
            %  @param as is the vector of polynomial cooefficients.
            %  @param t  is the domain of the polynomial.
            %  @returns ed, the range  of the polynomial.
            
            y = as(1)*ones('like', t);
            for ia = 2:length(as)
                y = y + as(ia)*t.^(ia-1);
            end
        end
    end  
    
	methods 
        function ed = estimateData(this, as)
            %% ESTIMATEDATA
            %  @param as is the vector of polynomial cooefficients.
            %  @returns ed, the estimated data isomorphic to this.dependentData.
            
            ed = this.polynomial(as, this.independentData);
        end
		  
 		function this = PolynomialsKernel(varargin)
 			%% POLYNOMIALSKERNEL
            
            this = this@mlbayesian.NullKernel(varargin{:});
 		end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

