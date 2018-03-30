classdef ConvolutionKernel < mlbayesian.NullKernel
	%% CONVOLUTIONKERNEL needs small data structures that will all fit in L1 caches. 

	%  $Revision$
 	%  was created 24-Dec-2017 17:52:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    methods (Static)
        function this = main(varargin)
            this = mlbayesian.ConvolutionKernel(varargin{:});
        end
        function y    = convolution(ps, t)
            %% CONVOLUTION
            %  @param ps is the vector [A T].
            %  @param t  is the domain of the convolution.
            %  @returns ed, the range  of the convolution.
            
            y = ps(1)*conv(exp(-t/ps(2)), exp(-t/ps(3)));
            y = y(1:length(t));
        end
    end  
    
	methods 
        function ed = estimateData(this, ps)
            %% ESTIMATEDATA
            %  @param ps is the vector of convolution cooefficients.
            %  @returns ed, the estimated data isomorphic to this.dependentData.
            
            ed = this.convolution(ps, this.independentData);
        end
		  
 		function this = ConvolutionKernel(varargin)
 			%% CONVOLUTIONKERNEL
            
            this = this@mlbayesian.NullKernel(varargin{:});		
 		end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

