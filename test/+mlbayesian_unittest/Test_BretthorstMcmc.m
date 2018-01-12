classdef Test_BretthorstMcmc < matlab.unittest.TestCase
	%% TEST_BRETTHORSTMCMC 

	%  Usage:  >> results = run(mlbayesian_unittest.Test_BretthorstMcmc)
 	%          >> result  = run(mlbayesian_unittest.Test_BretthorstMcmc, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 24-Dec-2017 15:29:01 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/test/+mlbayesian_unittest.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
        home = '/Users/jjlee/Tmp'
        pwd0
 		registry
 		testObj
 	end

	methods (Test)
		function test_linear(this)
            import mlbayesian.*;
            mdl = PolynomialsModel;
            mdl.as  = [1 1 0 0];
            mdl.sas = [1 1 0 0];
            mdl.fixed = logical([0 0 1 1]);
            mdl.fixedValue = [0 0 0 0];
            mdl = mdl.doConstructGenerative; % generate dependentData using as, sas, fixed, fixedValue listed above
            this.testObj.model = mdl;
            this.testObj = this.testObj.estimateParameters;
            this.testObj.diagnose;
            this.verifyTrue( this.testObj.isfinished);
            this.verifyEqual(this.testObj.model.modelParameters, mdl.modelParameters, 'RelTol', 1e-7);
            this.verifyEqual(this.testObj.model.objectiveFunc, 0, 'AbsTol', 1e-14);
            this.verifyEqual(this.testObj.model.fqfilename, fullfile(this.home, 'mlbayesian_PolynomialsModel.mat'));
 		end
		function test_polynomialsValidative(this)
            % See also:  Test_McmcCellular, test_polynomials
            
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            import mlbayesian.*;
            mdl = PolynomialsModel('independentData', t, 'dependentData', y);
            mdl.as  = [1 1 1 0]; % differing initial conditions; don't run do ConstructGenerative; compare to Test_McmcCellular
            mdl.sas = [1 1 1 0];
            mdl.asMin = [0 0 0 0];
            mdl.asMax = [10 10 10 0];
            mdl.fixed = logical([0 0 0 1]);
            mdl.fixedValue = [1 1 1 0];
            this.testObj.model = mdl;
            this.testObj = this.testObj.estimateParameters;
            this.testObj.diagnose;
            this.verifyTrue( this.testObj.isfinished);
            this.verifyEqual(this.testObj.model.modelParameters, [1 2 3 0], 'RelTol', 1e-4);
            this.verifyEqual(this.testObj.model.objectiveFunc, 0, 'AbsTol', 1e-5);
            this.verifyEqual(this.testObj.model.fqfilename, fullfile(this.home, 'mlbayesian_PolynomialsModel.mat'));
        end
		function test_polynomials(this)
            import mlbayesian.*;
            mdl = PolynomialsModel;
            mdl.nBeta = 100;
            mdl.nAnneal = 40;
            this.testObj.model = mdl;
            this.testObj = this.testObj.estimateParameters;
            this.testObj.diagnose;
            this.verifyTrue( this.testObj.isfinished);
            this.verifyEqual(this.testObj.model.modelParameters, mdl.modelParameters, 'RelTol', 5e-4);
            this.verifyEqual(this.testObj.model.objectiveFunc, 0, 'AbsTol', 1e-5);
            this.verifyEqual(this.testObj.model.fqfilename, fullfile(this.home, 'mlbayesian_PolynomialsModel.mat'));
        end
		function test_convolution(this)
            import mlbayesian.*;
            mdl = ConvolutionModel;
            mdl.nBeta = 100;
            mdl.nAnneal = 40;
            this.testObj.model = mdl;
            this.testObj = this.testObj.estimateParameters;
            this.testObj.diagnose;
            this.verifyTrue( this.testObj.isfinished);
            this.verifyEqual(this.testObj.model.modelParameters, mdl.modelParameters, 'RelTol', 5e-4);
            this.verifyEqual(this.testObj.model.objectiveFunc, 0, 'AbsTol', 1e-9);
            this.verifyEqual(this.testObj.model.fqfilename, fullfile(this.home, 'mlbayesian_ConvolutionModel.mat'));
 		end
	end

 	methods (TestClassSetup)
		function setupBretthorstMcmc(this)
 			import mlbayesian.*;
            this.pwd0 = pushd(this.home);
 			this.testObj_ = BretthorstMcmc;
            this.addTeardown(@this.cleanFiles);
 		end
	end

 	methods (TestMethodSetup)
		function setupBretthorstMcmcTest(this)
 			this.testObj = this.testObj_;
 		end
    end 
    
    %% PRIVATE

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
            popd(this.pwd0);
            close all
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

