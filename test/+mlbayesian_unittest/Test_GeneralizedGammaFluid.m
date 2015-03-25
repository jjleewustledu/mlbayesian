classdef Test_GeneralizedGammaFluid < matlab.unittest.TestCase 
	%% TEST_GENERALIZEDGAMMAFLUID  

	%  Usage:  >> results = run(mlbayesian_unittest.Test_GeneralizedGammaFluid)
 	%          >> result  = run(mlbayesian_unittest.Test_GeneralizedGammaFluid, 'test_dt')
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
        a   = 10
        d   = 3
        da  = 0
        dp  = 0
        fss = 0.1
        p   = 1
        q0  = 1
        t0  = 10
        times = 0:99
 	end 

	methods (Test) 
 		function test_plotAs(this)
            figure
            hold on
            a_ = [2 4 8 16 32 64];
            for idx = 1:length(a_)
                plot(this.times, ...
                     mlbayesian.GeneralizedGammaFluid.simulateQ(a_(idx), this.d, this.da, this.dp, this.fss, this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->var, d->%g, da->%g, dp->%g, fss->%g, p->%g, q0->%g, t0->%g', ...
                                                                         this.d, this.da, this.dp, this.fss, this.p, this.q0, this.t0));
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
                     mlbayesian.GeneralizedGammaFluid.simulateQ(this.a, d_(idx), this.da, this.dp, this.fss, this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->var, da->%g, dp->%g, fss->%g, p->%g, q0->%g, t0->%g', ...
                                                                this.a,          this.da, this.dp, this.fss, this.p, this.q0, this.t0));
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
                     mlbayesian.GeneralizedGammaFluid.simulateQ(this.a, this.d, this.da, this.dp, fss_(idx), this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, da->%g, dp->%g, fss->var, p->%g, q0->%g, t0->%g', ...
                                                                this.a, this.d, this.da, this.dp,            this.p, this.q0, this.t0));
            legend(cellfun(@(x) sprintf('fss = %g', x), num2cell(fss_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotPs(this)
            figure
            hold on
            p_ = [0.5 0.75 1 1.5 2 3 4 8];
            for idx = 1:length(p_)
                plot(this.times, ...
                     mlbayesian.GeneralizedGammaFluid.simulateQ(this.a, this.d, this.da, this.dp, this.fss, p_(idx), this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, da->%g, dp->%g, fss->%g, p->var, q0->%g, t0->%g', ...
                                                                this.a, this.d, this.da, this.dp, this.fss,          this.q0, this.t0));
            legend(cellfun(@(x) sprintf('p = %g', x), num2cell(p_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotQ0s(this)
            figure
            hold on
            q0_ = [0.5 1 2 4 8];
            for idx = 1:length(q0_)
                plot(this.times, ...
                     mlbayesian.GeneralizedGammaFluid.simulateQ(this.a, this.d, this.da, this.dp, this.fss, this.p, q0_(idx), this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, da->%g, dp->%g, fss->%g, p->%g, q0->var, t0->%g', ...
                                                                this.a, this.d, this.da, this.dp, this.fss, this.p,           this.t0));
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
                     mlbayesian.GeneralizedGammaFluid.simulateQ(this.a, this.d, this.da, this.dp, this.fss, this.p, this.q0, t0_(idx), this.times));
            end
            title(sprintf('a->%g, d->%g, da->%g, dp->%g, fss->%g, p->%g, q0->%g, t0->var', ...
                                                                this.a, this.d, this.da, this.dp, this.fss, this.p, this.q0));
            legend(cellfun(@(x) sprintf('t0 = %g', x), num2cell(t0_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
        end 
        function test_simulateMcmc(this)
            mlbayesian.GeneralizedGammaFluid.simulateMcmc(this.a, this.d, this.da, this.dp, this.fss, this.p, this.q0, this.t0, this.times);
        end
 	end 

 	methods (TestClassSetup) 
 		function setupGeneralizedGammaFluid(this) 
            import mlbayesian.*;
            this.testQ   = GeneralizedGammaFluid.simulateQ(this.a, this.d, this.da, this.dp, this.fss, this.p, this.q0, this.t0, this.times);
 			this.testObj = GeneralizedGammaFluid(this.times, this.testQ);
 		end 
 	end 

 	methods (TestClassTeardown) 
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

