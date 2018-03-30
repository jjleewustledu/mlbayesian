classdef Test_GeneralizedGammaStrategy < matlab.unittest.TestCase
	%% TEST_GENERALIZEDGAMMASTRATEGY 

	%  Usage:  >> results = run(mlbayesian_unittest.Test_GeneralizedGammaStrategy)
 	%          >> result  = run(mlbayesian_unittest.Test_GeneralizedGammaStrategy, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 27-Jun-2017 22:11:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/test/+mlbayesian_unittest.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
 		registry
 		testObj
 	end

	methods (Test)
        function test_godo(this)
            this.testObj = mlbayesian.GeneralizedGammaStrategy.godo;
        end
		function test_afun(this)
 			import mlbayesian.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_simulateItsMcmc(this)
            this.testObj.simulateItsMcmc;
        end
        function test_doItsBayes(this)
            this.testObj.doItsBayes;
        end
        function test_plotParVars(this)
            this.testObj.plotParVars( ...
                'a', [1 10 100]);
            this.testObj.plotParVars( ...
                'b', [1 10 100]);
            this.testObj.plotParVars( ...
                'p', [1 10 100]);
        end
	end

 	methods (TestClassSetup)
		function setupGeneralizedGammaStrategy(this)
 			import mlbayesian.*;
 			this.testObj_ = GeneralizedGammaStrategy;
 		end
	end

 	methods (TestMethodSetup)
		function setupGeneralizedGammaStrategyTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

