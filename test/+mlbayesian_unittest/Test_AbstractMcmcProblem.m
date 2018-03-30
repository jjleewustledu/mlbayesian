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
        unittest_home = fullfile(getenv('HOME'), 'MATLAB-Drive/mlbayesian/test/+mlbayesian_unittest')
        testProblem = 'kinetics4'
    end
    
	methods (Test)
 		function test_runMcmc(this)  
            disp(datestr(now))			
            switch (this.testProblem)
                case 'polynomial'
                    [y1,y] = this.polynomialCase;
                    this.verifyEqual(y1, y, 'RelTol', 0.001);
                case 'kinetics4'                    
                    [k1,k] = this.kinetics4Case;
                    this.verifyEqual(k1, k, 'RelTol', 0.05);
                case 'kinetics3'                    
                    [k1,k] = this.kinetics3Case;
                    this.verifyEqual(k1, k, 'RelTol', 0.05);
                otherwise
            end            
        end 
    end 
    
    methods 
        function [y1, y] = polynomialCase(~)
            t = 0:20;
            y = 1 + 2*t + 3*t.^2;
            
            pmp = mlbayesian_unittest.PolynomialMcmcProblem(t, y);
            pmp = pmp.estimateParameters;
            y1  = pmp.estimateData;
        end
        function [k1, k] = kineticsCase(this)
            k = [0.26057 0.00047967 0.034772 0.0025173 0.00042528 4.3924]; % k04 k12 k21 k32 k43 t0
            cd(this.unittest_home);
            load('fourCompSimTimeInterpolants');
            load('fourCompSimQ');  
            load('fourCompSimCa');
            
            kmp = mlbayesian_unittest.Kinetics4McmcProblem( ...
                  fourCompSimTimeInterpolants, fourCompSimQ, fourCompSimCa);
            kmp = kmp.estimateParameters;
            k1  = [kmp.finalParams('k04'), kmp.finalParams('k12'), kmp.finalParams('k21'), ...
                   kmp.finalParams('k32'), kmp.finalParams('k43'), kmp.finalParams('t0')];
            figure; plot([fourCompSimQ' kmp.estimateData'])
        end 
        function [k1, k] = kinetics4Case(this)             
            cd(this.unittest_home);
            import mlpet.*;
            dta = DTA.load('p5661g.dta');
            tsc = TSC.import('p5661wb.tsc');
            len = min(length(dta.timeInterpolants), length(tsc.timeInterpolants));
            timeInterp = tsc.timeInterpolants(1:len);
            Ca = dta.countInterpolants(1:len);
            Q = tsc.activityInterpolants(1:len);            
            figure; plot(timeInterp, Ca, timeInterp, Q)
            
            kmp = mlbayesian_unittest.Kinetics4McmcProblem( ...
                  timeInterp, Q, Ca);
            kmp = kmp.estimateParameters;
            k   = [kmp.finalParams('k04'), kmp.finalParams('k12'), kmp.finalParams('k21'), ...
                   kmp.finalParams('k32'), kmp.finalParams('k43'), kmp.finalParams('t0')];   
            figure; plot(timeInterp, Ca, timeInterp, kmp.estimateData)         
            save('Test_AbstractMcmcProblem_kinetics4Case_k.mat', 'k');
            
            kmp1 = mlbayesian_unittest.Kinetics4McmcProblem( ...
                   timeInterp, kmp.estimateData, Ca);
            kmp1 = kmp1.estimateParameters;
            k1   = [kmp1.finalParams('k04'), kmp1.finalParams('k12'), kmp1.finalParams('k21'), ...
                    kmp1.finalParams('k32'), kmp1.finalParams('k43'), kmp1.finalParams('t0')];
            figure; plot(timeInterp, Ca, timeInterp, kmp.estimateData, timeInterp, kmp1.estimateData)
        end 
        function [k1, k] = kinetics3Case(this)             
            cd(this.unittest_home);
            import mlpet.*;
            dta = DTA.load('p5661g.dta');
            tsc = TSC.import('p5661wb.tsc');
            len = min(dta.scanDuration, tsc.scanDuration);
            timeInterp = tsc.timeInterpolants(1:len);
            Ca = dta.countInterpolants(1:len);
            Q = tsc.activityInterpolants(1:len);            
            figure; plot(timeInterp, Ca, timeInterp, Q)
            
            kmp = mlbayesian_unittest.Kinetics3McmcProblem( ...
                  timeInterp, Q, Ca);
            kmp = kmp.estimateParameters;
            k   = [kmp.finalParams('k04'), kmp.finalParams('k12'), kmp.finalParams('k21'), ...
                   kmp.finalParams('k32'),                         kmp.finalParams('t0')];   
            figure; plot(timeInterp, Ca, timeInterp, kmp.estimateData)         
            save('Test_AbstractMcmcProblem_kinetics3Case_k.mat', 'k');
            
            kmp1 = mlbayesian_unittest.Kinetics3McmcProblem( ...
                   timeInterp, kmp.estimateData, Ca);
            kmp1 = kmp1.estimateParameters;
            k1   = [kmp1.finalParams('k04'), kmp1.finalParams('k12'), kmp1.finalParams('k21'), ...
                    kmp1.finalParams('k32'),                          kmp1.finalParams('t0')];
            figure; plot(timeInterp, Ca, timeInterp, kmp.estimateData, timeInterp, kmp1.estimateData)
        end 
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

