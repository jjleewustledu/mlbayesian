classdef BestGammaFluid < mlbayesian.AbstractMcmcProblem
	%% BESTGAMMAFLUID fits a generalized gamma & exponential recovery to time-counts data.
    %  http://en.wikipedia.org/wiki/Generalized_gamma_distribution
    %  N.B.  f(tau; a,d,p) = \Gamma^{-1}(d/p) (p/a^d) tau^(d-1) exp(-(tau/a)^p) with a > 0, d > 0, p > 0, t - t0 > 0.

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 

	properties
        showPlots = true	 
        baseTitle = 'BestGammaFluid'
        xLabel    = 'times/s'
        yLabel    = 'wellcounts'
    end 
    
    methods (Static)  
        function Q    = simulateQ(a, d, da, dp, fss, p, q0, t0, times)
            idx_t0 = floor(t0) + 1;            
            cnorm  = q0 * ((p/a^d)/gamma(d/p));
            exp1   = abs(exp(-(times/a).^p));
            exp2   = abs(exp(-(times/(a+da)).^(p+dp)));
            Q0     = cnorm * times.^(d-1) .* exp1;
            Q0     = Q0 + max(Q0) * fss * (1 - exp2); 
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
            assert(all(isreal(Q)), 'BestGammaFluid.simulateQ.Q was complex');
            assert(~any(isnan(Q)), 'BestGammaFluid.simulateQ.Q was NaN: %s', num2str(Q));
        end
        function this = simulateMcmc(a, d, da, dp, fss, p, q0, t0, times)
            
            import mlbayesian.*;            
            Q    = BestGammaFluid.simulateQ(a, d, da, dp, fss, p, q0, t0, times);
            this = BestGammaFluid(times, Q);
            this = this.estimateParameters %#ok<NOPRT>
            
            figure;
            plot(times, this.estimateData, times, Q, 'o');
            legend('Bayesian estimate', 'simulated');
            title(sprintf('simulateMcmc expected:  a->%g, d->%g, da->%g, dp->%g, fss->%g, p->%g, q0->%g, t0->%g, max(t)->%g', ...
                  a, d, da, dp, fss, p, q0, t0, max(times)));
            xlabel('time/s');
            ylabel('generalized gamma/arbitrary');
        end
    end

	methods 		  
 		function this = BestGammaFluid(varargin) 
 			%% BESTGAMMAFLUID 
 			%  Usage:  this = BestGammaFluid(times, well_counts) 

 			this = this@mlbayesian.AbstractMcmcProblem(varargin{:}); 
 		end 
        function this = estimateParameters(this)
            %% ESTIMATEPARAMETERS a, d, da, dp, p, q0, t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('a')   = struct('fixed', 0, 'min', eps, 'mean',  8.5,   'max', this.timeFinal/2);
            map('d')   = struct('fixed', 0, 'min', eps, 'mean',  5.4,   'max',    8);
            map('fss') = struct('fixed', 0, 'min',   0, 'mean',  0.46,  'max',    0.66);
            map('p')   = struct('fixed', 0, 'min', eps, 'mean',  1.1,   'max',    2); 
            map('q0')  = struct('fixed', 0, 'min', eps, 'mean',  3.7e6, 'max',   50*max(this.dependentData));
            map('t0')  = struct('fixed', 0, 'min',   0, 'mean', 31,     'max', this.timeFinal/2);            
            
            map('da')  = struct('fixed', 1, 'min', -7.5, 'mean', 0, 'max', this.timeFinal/2);
            map('dp')  = struct('fixed', 1, 'min', -1,   'mean', 0, 'max', 1);

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('a'),  this.finalParams('d'), this.finalParams('da'), this.finalParams('dp'),...
                this.finalParams('fss'),this.finalParams('p'), this.finalParams('q0'), this.finalParams('t0'));
        end
        function Q    = estimateDataFast(this, a, d, da, dp, fss, p, q0, t0)  
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            cnorm  = q0 * ((p/a^d)/gamma(d/p));
            exp1   = abs(exp(-(times/a).^p));
            exp2   = abs(exp(-(times/(a+da)).^(p+dp)));
            Q0     = cnorm * times.^(d-1) .* exp1;
            Q0     = Q0 + max(Q0) * fss * (1 - exp2);    
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
            assert(all(isreal(Q)), 'BestGammaFluid.estimateDataFast.Q was complex');
            assert(~any(isnan(Q)), 'BestGammaFluid.estimateDataFast.Q was NaN: %s', num2str(Q));
        end  
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end
