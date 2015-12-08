classdef Test_McmcCellular < matlab.unittest.TestCase
	%% TEST_MCMCCELLULAR  
    %
    %  Usage:  >> results = run(mlbayesian_unittest.Test_McmcCellular)
	%          >> result  = run(mlbayesian_unittest.Test_McmcCellular, 'test_dt')
	%  See also:  file:///Applications/Developer/MATLAB_R2014a.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 23-Nov-2015 17:37:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	 
	methods (Test)
 		function test_polynomial(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y});
            pmc = pmc.estimateParameters;
            y1  = pmc.estimateData;
            
            this.verifyEqual(y1{1}, y, 'RelTol', 0.001);
 		end 
    end 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

