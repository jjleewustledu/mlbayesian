classdef Test_MCMC < matlab.unittest.TestCase
	%% TEST_MCMC  
    %
    %  Usage:  >> results = run(mlbayesian_unittest.Test_MCMC)
	%          >> result  = run(mlbayesian_unittest.Test_MCMC, 'test_dt')
	%  See also:  file:///Applications/Developer/MATLAB_R2014a.app/help/matlab/matlab-unit-test-framework.html
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
 	 
    properties 
        testProblem = 'polynomial';
    end
    
	methods (Test)
 		function test_runMcmc(this)  			
            switch (this.testProblem)
                case 'polynomial'
                    [y1,y] = this.polynomialCase;
                    this.verifyEqual(y1, y, 'RelTol', 0.001);
                otherwise
            end            
 		end 
    end 

    %% PROTECTED
    
    methods (Access = 'protected')
        function [y1, y] = polynomialCase(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pp = mlbayesian_unittest.PolynomialProblem(t, y);
            pp = pp.estimateParameters;
            y1 = pp.estimateData;
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

