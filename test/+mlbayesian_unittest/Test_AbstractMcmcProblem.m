classdef Test_AbstractMcmcProblem < matlab.unittest.TestCase 
	%% TEST_ABSTRACTMCMCPROBLEM  

	%  Usage:  >> results = run(mlbayesian_unittest.Test_AbstractMcmcProblem)
 	%          >> result  = run(mlbayesian_unittest.Test_AbstractMcmcProblem, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 

    properties 
        unittest_home = '/Users/jjlee/Local/src/mlcvl/mlbayesian/test/+mlbayesian_unittest'
        testProblem = 'kinetics'
    end
    
	methods (Test)
 		function test_runMcmc(this)  			
            switch (this.testProblem)
                case 'polynomial'
                    [y1,y] = this.polynomialCase;
                    this.verifyEqual(y1, y, 'RelTol', 0.001);
                case 'kinetics'                    
                    [k1,k] = this.kineticsCase;
                    this.verifyEqual(k1, k, 'RelTol', 0.05);
                otherwise
            end            
        end 
    end 

    %% PROTECTED
    
    methods (Access = 'protected')
        function [y1, y] = polynomialCase(~)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmp = mlbayesian_unittest.PolynomialMcmcProblem(t, y);
            pmp = pmp.estimateParameters;
            y1  = pmp.estimateData;
        end
        function [k1, k] = kineticsCase(this)
            k = [0.26057 0.00047967 0.034772 0.0025173 0.00042528 0]; % k04 k12 k21 k32 k43 t0
            cd(this.unittest_home);
            load('fourCompartmentsSimulatorTimeInterpolants');
            load('fourCompartmentsSimulatorQ');  
            load('fourCompartmentsSimulatorCa');
            
            kmp = mlbayesian_unittest.KineticsMcmcProblem( ...
                  fourCompartmentsSimulatorTimeInterpolants, fourCompartmentsSimulatorQ, fourCompartmentsSimulatorCa);
            kmp = kmp.estimateParameters;
            k1  = [kmp.finalParams('k04'), kmp.finalParams('k12'), kmp.finalParams('k21'), ...
                   kmp.finalParams('k32'), kmp.finalParams('k43'), kmp.finalParams('t0')];
            figure; plot([fourCompartmentsSimulatorQ' kmp.estimateData'])
        end 
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

