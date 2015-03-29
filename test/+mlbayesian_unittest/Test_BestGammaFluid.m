classdef Test_BestGammaFluid < matlab.unittest.TestCase 
	%% TEST_BESTGAMMAFLUID  

	%  Usage:  >> results = run(mlbayesian_unittest.Test_BestGammaFluid)
 	%          >> result  = run(mlbayesian_unittest.Test_BestGammaFluid, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 

	properties 
        testQ
 		testObj 
        a   =   8.5
        d   =   5.4
        fss =   0.46
        p   =  1.1
        q0  =  3.7e6
        t0  = 31
        times = 0:119
        testFolder = '/Users/jjlee/Local/src/mlcvl/mlbayesian/test/+mlbayesian_unittest'
    end 
    
    properties (Dependent)
        expectedBestFitParams
    end
    
    methods %% GET/SET
        function e = get.expectedBestFitParams(this)
            e = [this.a this.d this.fss this.p this.q0 this.t0]';
        end
    end

	methods (Test) 
 		function test_plotAs(this)
            figure
            hold on
            a_ = [1 2 4 8 16 32];
            for idx = 1:length(a_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.simulateQ(a_(idx), this.d, this.fss, this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->var, d->%g, fss->%g, p->%g, q0->%g, t0->%g', ...
                                                                  this.d, this.fss, this.p, this.q0, this.t0));
            legend(cellfun(@(x) sprintf('a = %g', x), num2cell(a_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotDs(this)
            figure
            hold on
            d_ = [1 2 3 4 5 6 7 8 9];
            for idx = 1:length(d_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.simulateQ(this.a, d_(idx), this.fss, this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->var, fss->%g, p->%g, q0->%g, t0->%g', ...
                                                         this.a,          this.fss, this.p, this.q0, this.t0));
            legend(cellfun(@(x) sprintf('d = %g', x), num2cell(d_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end  
 		function test_plotFSSs(this)
            figure
            hold on
            fss_ = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
            for idx = 1:length(fss_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.simulateQ(this.a, this.d, fss_(idx), this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, fss->var, p->%g, q0->%g, t0->%g', ...
                                                         this.a, this.d,            this.p, this.q0, this.t0));
            legend(cellfun(@(x) sprintf('fss = %g', x), num2cell(fss_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotPs(this)
            figure
            hold on
            p_ = [0.5 0.75 0.9 1 1.1 1.5 2 3];
            for idx = 1:length(p_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.simulateQ(this.a, this.d, this.fss, p_(idx), this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, fss->%g, p->var, q0->%g, t0->%g', ...
                                                         this.a, this.d, this.fss,          this.q0, this.t0));
            legend(cellfun(@(x) sprintf('p = %g', x), num2cell(p_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotQ0s(this)
            figure
            hold on
            q0_ = [0.5 1 2 4 8] * 1e7;
            for idx = 1:length(q0_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.simulateQ(this.a, this.d, this.fss, this.p, q0_(idx), this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, fss->%g, p->%g, q0->var, t0->%g', ...
                                                         this.a, this.d, this.fss, this.p,           this.t0));
            legend(cellfun(@(x) sprintf('q0 = %g', x), num2cell(q0_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotT0s(this)
            figure
            hold on
            t0_ = [0 4 8 16 32 64];
            for idx = 1:length(t0_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.simulateQ(this.a, this.d, this.fss, this.p, this.q0, t0_(idx), this.times));
            end
            title(sprintf('a->%g, d->%g, fss->%g, p->%g, q0->%g, t0->var', ...
                                                         this.a, this.d, this.fss, this.p, this.q0));
            legend(cellfun(@(x) sprintf('t0 = %g', x), num2cell(t0_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
        end 
        function test_simulateMcmc(this)
            this.testObj = mlbayesian.BestGammaFluid.simulateMcmc(this.a, this.d, this.fss, this.p, this.q0, this.t0, this.times);
            this.assertEqual(this.testObj.bestFitParams, this.expectedBestFitParams, 'RelTol', 0.05);
            bestGammaFluid = this.testObj; %#ok<NASGU>
            save(fullfile(this.testFolder, 'bestGammaFluid.mat'), 'bestGammaFluid');
        end
 	end 

 	methods (TestClassSetup) 
 		function setupBestGammaFluid(this) 
            import mlbayesian.*;
            this.testQ   = BestGammaFluid.simulateQ(this.a, this.d, this.fss, this.p, this.q0, this.t0, this.times);
 			this.testObj = BestGammaFluid(this.times, this.testQ);
 		end 
 	end 

 	methods (TestClassTeardown) 
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

