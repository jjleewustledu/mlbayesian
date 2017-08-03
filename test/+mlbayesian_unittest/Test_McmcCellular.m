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
            
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'mcmcParameters', this.mcmcParameters, ...
                'NyquistFreqFactor', 2);
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 0.0005);
        end
 		function test_dt(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'mcmcParameters', this.mcmcParameters, ...
                'NyquistFreqFactor', 8);
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 5e-8);
            this.verifyEqual(pmc.dt, 1/8, 'AbsTol', eps);
            this.verifyEqual(length(pmc.independentData{1}), 21);
            this.verifyEqual(length(pmc.dependentData{1}),   21);
            this.verifyEqual(length(pmc.independentDataInterp{1}), 161);
            this.verifyEqual(length(pmc.dependentDataInterp{1}),   161);
        end         
 		function test_hist3(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'mcmcParameters', this.mcmcParameters, ...
                'NyquistFreqFactor', 2);
            pmc = pmc.estimateParameters;
        end 
        function test_polynomial_dProposals(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            mp = this.mcmcParameters;
            mp.nProposals = 400; % 20% increase in run time, rising plotLogProbabilityQC
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'mcmcParameters', mp, ...
                'NyquistFreqFactor', 2);
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 0.0005);
        end
        function test_polynomial_dPop(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            mp = this.mcmcParameters;
            mp.nPop = 200; % 380% increase in run time
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'mcmcParameters', mp, ...
                'NyquistFreqFactor', 2);
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 0.0005);
        end
        function test_polynomial_dBeta(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            mp = this.mcmcParameters;
            mp.nBeta = 200; % 358% incr run time; 2x wider histStdOfError; 10^{-4} decr parameter posterior widths; plateau of plotAnnealing
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'mcmcParameters', mp, ...
                'NyquistFreqFactor', 2);
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 5e-8);
        end
        function test_polynomial_dAnneal2(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'NyquistFreqFactor', 2);
            pmc = pmc.adjustN('nAnneal', 2); % 187% incr run time; 2x wider histStdOfError; 10^{-3} decr parameter posterior widths
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 5e-7);
        end
        function test_polynomial_dAnneal1_5(this)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmc = mlbayesian_unittest.PolynomialMcmcCellular({t}, {y}, ...
                'NyquistFreqFactor', 2);
            pmc = pmc.adjustN('nAnneal', 1.5); % 137% incr run time; 10^{-1} decr parameter posterior widths
            tic; pmc = pmc.estimateParameters; toc
            y1 = pmc.estimateData;            
            this.verifyEqual(y1{1}, y, 'RelTol', 5e-5);
        end
    end 
    
    methods 
        function mcmcp = mcmcParameters(this)
            m = containers.Map;
            m('c0') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            m('c1') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            m('c2') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10);
            mcmcp = mlbayesian.McmcParameters(m, 21);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

