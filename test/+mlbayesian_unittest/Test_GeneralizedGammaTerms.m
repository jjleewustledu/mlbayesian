classdef Test_GeneralizedGammaTerms < matlab.unittest.TestCase
	%% TEST_GENERALIZEDGAMMATERMS 

	%  Usage:  >> results = run(mlbayesian_unittest.Test_GeneralizedGammaTerms)
 	%          >> result  = run(mlbayesian_unittest.Test_GeneralizedGammaTerms, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 27-Jun-2017 22:38:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/test/+mlbayesian_unittest.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
 		registry
 		testObj
        t = 0:120
 	end

	methods (Test)
		function test_afun(this)
 			import mlbayesian.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_steadyState(this)
            figure;            
            hold on;
            plot(this.testObj.steadyState(1, 1/6,  0, this.t));
            plot(this.testObj.steadyState(1, 1/12, 0, this.t));
            plot(this.testObj.steadyState(1, 1/24, 0, this.t));
            plot(this.testObj.steadyState(1, 1/36, 0, this.t));
            plot(this.testObj.steadyState(1, 1/48, 0, this.t));
            plot(this.testObj.steadyState(1, 1/64, 0, this.t));
            xlabel('time/s');
            ylabel('steadyState')
            title('S = 1, t0 = 0');
            legend('k = 1/6', ' 1/12', '1/24', '1/36', '1/48', '1/64');
        end
        function test_gammaTerm(this)
            figure;            
            hold on;
            plot(this.testObj.gammaTerm(1/32, 1/36, 12, this.t));
            plot(this.testObj.gammaTerm(1/16, 1/36, 12, this.t));
            plot(this.testObj.gammaTerm(1/4, 1/36, 12, this.t));
            plot(this.testObj.gammaTerm(1/2, 1/36, 12, this.t));
            plot(this.testObj.gammaTerm(1,   1/36, 12, this.t));
            plot(this.testObj.gammaTerm(2,   1/36, 12, this.t));
            plot(this.testObj.gammaTerm(4,   1/36, 12, this.t));
            plot(this.testObj.gammaTerm(16,  1/36, 12, this.t));
            plot(this.testObj.gammaTerm(32,  1/36, 12, this.t));
            xlabel('time/s');
            ylabel('gammaTerm')
            title('b = 1/36; t0 = 12');
            legend('a = 1/32', '1/16', '1/4', '1/2', '1', '2', '4', '16', '32');
            
            figure;            
            hold on;
            plot(this.testObj.gammaTerm(1.5, 12,   12, this.t));
            plot(this.testObj.gammaTerm(1.5, 6,    12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1,    12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/6,  12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/12, 12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/24, 12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36, 12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/48, 12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/64, 12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/128, 12, this.t));
            xlabel('time/s');
            ylabel('gammaTerm')
            title('a = 1.5; t0 = 12');
            legend('b = 12', '6', '1', '1/6', '1/12', '1/24', '1/36', '1/48', '1/64', '1/128');
            
            figure;            
            hold on;
            plot(this.testObj.gammaTerm(1.5, 1/36, -64, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36, -48, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36, -36, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36, -24, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36, -12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36,   0, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36,  12, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36,  24, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36,  36, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36,  48, this.t));
            plot(this.testObj.gammaTerm(1.5, 1/36,  64, this.t));
            xlabel('time/s');
            ylabel('gammaTerm')
            title('a = 1.5; b = 1/6');
            legend('t0 = -64', '-48', '-36', '-24', '-12', '0', '12', '24', '36', '48', '64');
        end
        function test_gammaSeries(this)
            figure;            
            hold on;
            plot(this.testObj.gammaSeries(1.5, 1/12, 0, 1.5, 1/12, 6,  0.5, this.t));
            plot(this.testObj.gammaSeries(1.5, 1/12, 0, 1.5, 1/12, 12, 0.5, this.t));
            plot(this.testObj.gammaSeries(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, this.t));
            plot(this.testObj.gammaSeries(1.5, 1/12, 0, 1.5, 1/12, 36, 0.5, this.t));
            plot(this.testObj.gammaSeries(1.5, 1/12, 0, 1.5, 1/12, 48, 0.5, this.t));
            plot(this.testObj.gammaSeries(1.5, 1/12, 0, 1.5, 1/12, 64, 0.5, this.t));
            xlabel('time/s');
            ylabel('gammaSeries')
            title('a = 1.5, b = 1/12, t01 = 0, weight = 0.5');
            legend('t02 = 6', '12', '24', '36', '48', '64');
        end
        function test_gammaSeriesSteady(this)
            figure;            
            hold on;
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/6,  24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/12, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/36, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/48, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/64, 24, this.t));
            xlabel('time/s');
            ylabel('gammaSeriesSteady')
            title('a = 1.5, b = 1/12, t01 = 0, t02 = 24, weight = 0.5, S = 0.005, t0 = 24');
            legend('k = 1/6', 'k = 1/12', 'k = 1/24', 'k = 1/36', 'k = 1/48', 'k = 1/64');
            
            figure;            
            hold on;
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.002, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.004, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.006, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.008, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.01,  1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.012, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.014, 1/24, 24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.016, 1/24, 24, this.t));
            xlabel('time/s');
            ylabel('gammaSeriesSteady')
            title('a = 1.5, b = 1/12, t01 = 0, t02 = 24, weight = 0.5, k = 1/24, t0 = 24');
            legend('S = 0.002', '0.004', '0.006', '0.008', '0.01', '0.012', '0.014', '0.016');
            
            figure;            
            hold on;
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24, -32, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24, -24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24, -16, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24,  -8, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24,   0, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24,   8, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24,  16, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24,  24, this.t));
            plot(this.testObj.gammaSeriesSteady(1.5, 1/12, 0, 1.5, 1/12, 24, 0.5, 0.005, 1/24,  32, this.t));
            xlabel('time/s');
            ylabel('gammaSeriesSteady')
            title('a = 1.5, b = 1/12, t01 = 0, t02 = 24, weight = 0.5, S = 0.005, k = 1/24');
            legend('t0 = -32', '-24', '-16', '-8', '0', '8', '16', '24', '32');
        end
        function test_gammaStretch(this)
            figure;            
            hold on;
            plot(this.testObj.gammaStretch(1.5, 1/6, 1/4, 0, this.t));
            plot(this.testObj.gammaStretch(1.5, 1/6, 1/3, 0, this.t));
            plot(this.testObj.gammaStretch(1.5, 1/6, 1/2, 0, this.t));
            plot(this.testObj.gammaStretch(1.5, 1/6, 1,   0, this.t));
            plot(this.testObj.gammaStretch(1.5, 1/6, 2,   0, this.t));
            plot(this.testObj.gammaStretch(1.5, 1/6, 3,   0, this.t));
            plot(this.testObj.gammaStretch(1.5, 1/6, 4,   0, this.t));
            xlabel('time/s');
            ylabel('gammaStretch')
            title('a = 1.5; b = 1/6, t0 = 0');
            legend('p = 1/4', '1/3', '1/2', '1', '2', '3', '4');
        end
        function test_gammaStretchSeries(this)
            figure;            
            hold on;
            plot(this.testObj.gammaStretchSeries(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 6,  0.5, this.t));
            plot(this.testObj.gammaStretchSeries(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 12, 0.5, this.t));
            plot(this.testObj.gammaStretchSeries(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, this.t));
            plot(this.testObj.gammaStretchSeries(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 36, 0.5, this.t));
            plot(this.testObj.gammaStretchSeries(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 48, 0.5, this.t));
            plot(this.testObj.gammaStretchSeries(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 64, 0.5, this.t));
            xlabel('time/s');
            ylabel('gammaStretchSeries')
            title('a = 1.5, b = 1/6, p = 3/2, t01 = 0, weight = 0.5');
            legend('t02 = 6', '12', '24', '36', '48', '64');
        end
        function test_gammaStretchSeriesSteady(this)
            figure;            
            hold on;
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24, -32, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24, -24, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24, -16, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24,  -8, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24,   0, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24,   8, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24,  16, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24,  24, this.t));
            plot(this.testObj.gammaStretchSeriesSteady(1.5, 1/6, 3/2, 0, 1.5, 1/6, 3/2, 24, 0.5, 0.005, 1/24,  32, this.t));
            xlabel('time/s');
            ylabel('gammaStretchSeriesSteady')
            title('a = 1.5, b = 1/6, p = 3/2, t01 = 0, t02 = 24, weight = 0.5, S = 0.005, k = 1/24');
            legend('t0 = -32', '-24', '-16', '-8', '0', '8', '16', '24', '32');
        end
	end

 	methods (TestClassSetup)
		function setupGeneralizedGammaTerms(this)
 			import mlbayesian.*;
 			this.testObj_ = GeneralizedGammaTerms;
 		end
	end

 	methods (TestMethodSetup)
		function setupGeneralizedGammaTermsTest(this)
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

