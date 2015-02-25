classdef Test_AbstractBayesianProblem < matlab.unittest.TestCase 
	%% TEST_ABSTRACTBAYESIANPROBLEM  

	%  Usage:  >> results = run(mlbayesian_unittest.Test_AbstractBayesianProblem)
 	%          >> result  = run(mlbayesian_unittest.Test_AbstractBayesianProblem, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
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
        function [y1, y] = polynomialCase(~)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pbp = mlbayesian_unittest.PolynomialBayesianProblem(t, y);
            pbp = pbp.estimateParameters;
            y1  = pbp.estimateData;
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

