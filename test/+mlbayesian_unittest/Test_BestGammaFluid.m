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
        testFluidQ
 		testObj 
        times = 0:119
        testFolder = '/Users/jjlee/Local/src/mlcvl/mlbayesian/test/+mlbayesian_unittest'
    end 
    
    properties (Dependent)        
        a
        d
        e
        p
        q0
        t0
    end
    
    methods %% GET
        function x = get.a(this)
            x = this.testObj.a;
        end
        function x = get.d(this)
            x = this.testObj.d;
        end
        function x = get.e(this)
            x = this.testObj.e;
        end
        function x = get.p(this)
            x = this.testObj.p;
        end
        function x = get.q0(this)
            x = this.testObj.q0;
        end
        function x = get.t0(this)
            x = this.testObj.t0;
        end
    end

	methods (Test) 
 		function test_plotAs(this)
            figure
            hold on
            a_ = [1 2 4 8 16 32];
            for idx = 1:length(a_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.fluidQ(a_(idx), this.d, this.e, this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->var, d->%g, e->%g, p->%g, q0->%g, t0->%g', ...
                                                                  this.d, this.e, this.p, this.q0, this.t0));
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
                     mlbayesian.BestGammaFluid.fluidQ(this.a, d_(idx), this.e, this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->var, e->%g, p->%g, q0->%g, t0->%g', ...
                                                         this.a,          this.e, this.p, this.q0, this.t0));
            legend(cellfun(@(x) sprintf('d = %g', x), num2cell(d_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end  
 		function test_plotEs(this)
            figure
            hold on
            e_ = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
            for idx = 1:length(e_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.fluidQ(this.a, this.d, e_(idx), this.p, this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, e->var, p->%g, q0->%g, t0->%g', ...
                                                         this.a, this.d,            this.p, this.q0, this.t0));
            legend(cellfun(@(x) sprintf('e = %g', x), num2cell(e_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
 		end 
 		function test_plotPs(this)
            figure
            hold on
            p_ = [0.5 0.75 0.9 1 1.1 1.5 2 3];
            for idx = 1:length(p_)
                plot(this.times, ...
                     mlbayesian.BestGammaFluid.fluidQ(this.a, this.d, this.e, p_(idx), this.q0, this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, e->%g, p->var, q0->%g, t0->%g', ...
                                                         this.a, this.d, this.e,          this.q0, this.t0));
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
                     mlbayesian.BestGammaFluid.fluidQ(this.a, this.d, this.e, this.p, q0_(idx), this.t0, this.times));
            end
            title(sprintf('a->%g, d->%g, e->%g, p->%g, q0->var, t0->%g', ...
                                                         this.a, this.d, this.e, this.p,           this.t0));
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
                     mlbayesian.BestGammaFluid.fluidQ(this.a, this.d, this.e, this.p, this.q0, t0_(idx), this.times));
            end
            title(sprintf('a->%g, d->%g, e->%g, p->%g, q0->%g, t0->var', ...
                                                         this.a, this.d, this.e, this.p, this.q0));
            legend(cellfun(@(x) sprintf('t0 = %g', x), num2cell(t0_), 'UniformOutput', false));
            xlabel('time/s');
            ylabel('arbitrary');
        end 
        
        function test_simulateMcmc(this)
            import mlbayesian.*;
            this.testObj = BestGammaFluid;            
            aMap = containers.Map; 
            fL = 0.85; fH = 1.15;  
            aMap('a')  = struct('fixed', 0, 'min', fL*this.a,  'mean', this.a,  'max', fH*this.a); 
            aMap('d')  = struct('fixed', 0, 'min', fL*this.d,  'mean', this.d,  'max', fH*this.d);
            aMap('e')  = struct('fixed', 0, 'min', fL*this.e,  'mean', this.e,  'max', fH*this.e);
            aMap('p')  = struct('fixed', 0, 'min', fL*this.p,  'mean', this.p,  'max', fH*this.p);
            aMap('q0') = struct('fixed', 0, 'min', fL*this.q0, 'mean', this.q0, 'max', fH*this.q0); 
            aMap('t0') = struct('fixed', 0, 'min', fL*this.t0, 'mean', this.t0, 'max', fH*this.t0); 
            
            o = BestGammaFluid.simulateMcmc(this.a, this.d, this.e, this.p, this.q0, this.t0, this.times, aMap);
            this.assertEqual(o.bestFitParams, o.expectedBestFitParams, 'RelTol', 0.05);
        end
        function test_GammaFluid(this)            
            this.testObj = mlbayesian.BestGammaFluid.runGammaFluid(this.times, this.testFluidQ);
            o = this.testObj;
            
            figure;
            plot(o.times, o.estimateData, this.times, this.testFluidQ, 'o');
            legend('Bayesian estimate', 'simulated data');
            title(sprintf('Laif1:  a %g, d %g, e %g, p %g, q0 %g, t0 %g', ...
                  o.a, o.d, o.e, o.p, o.q0, o.t0));
            xlabel('time/s');
            ylabel('arbitrary');
            
            this.assertEqual(o.bestFitParams, o.expectedBestFitParams, 'RelTol', 0.05);
        end
 	end 

 	methods (TestClassSetup) 
 		function this = setupBestGammaFluid(this)
            this.testObj = mlbayesian.BestGammaFluid; 
 		end 
 	end 

 	methods (TestClassTeardown) 
    end 
    
    methods
        function this = Test_BestGammaFluid
             this = this@matlab.unittest.TestCase;
             import mlbayesian.*;
             f = BestGammaFluid;
             this.testFluidQ = BestGammaFluid.fluidQ(f.a, f.d, f.e, f.p, f.q0, f.t0, this.times);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

