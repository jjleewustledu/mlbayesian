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
        function Q    = simulateQ(a, d, fss, p, q0, t0, times)
            idx_t0 = floor(t0) + 1;            
            cnorm  = q0 * ((p/a^d)/gamma(d/p));
            exp1   = abs(exp(-(times/a).^p));
            Q0     = cnorm * times.^(d-1) .* exp1;
            Q0     = Q0 + max(Q0) * fss * (1 - exp1); 
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
            assert(all(isreal(Q)), 'BestGammaFluid.simulateQ.Q was complex');
            assert(~any(isnan(Q)), 'BestGammaFluid.simulateQ.Q was NaN: %s', num2str(Q));
        end
        function this = simulateMcmc(a, d, fss, p, q0, t0, times)
            
            import mlbayesian.*;            
            Q    = BestGammaFluid.simulateQ(a, d, fss, p, q0, t0, times);
            this = BestGammaFluid(times, Q);
            this = this.estimateParameters %#ok<NOPRT>
            
            figure;
            plot(times, this.estimateData, times, Q, 'o');
            legend('Bayesian estimate', 'simulated');
            title(sprintf('simulateMcmc expected:  a->%g, d->%g, fss->%g, p->%g, q0->%g, t0->%g, max(t)->%g', ...
                  a, d, fss, p, q0, t0, max(times)));
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
            %% ESTIMATEPARAMETERS a, d, fss, p, q0, t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('a')   = struct('fixed', 0, 'min',   2,   'mean',  8.5,   'max',   16);
            map('d')   = struct('fixed', 0, 'min',   4,   'mean',  5.4,   'max',    8);
            map('fss') = struct('fixed', 0, 'min',   0.3, 'mean',  0.46,  'max',    0.6);
            map('p')   = struct('fixed', 0, 'min',   0.5, 'mean',  1.1,   'max',    1.5); 
            map('q0')  = struct('fixed', 0, 'min', 1e5,   'mean',  3.7e6, 'max',   50*max(this.dependentData));
            map('t0')  = struct('fixed', 0, 'min',   0,   'mean', 31,     'max', this.timeFinal/2); 

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function this = estimateDcvParameters(this)
            %% ESTIMATEDCVPARAMETERS a, d, fss, p, q0, t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('a')   = struct('fixed', 0, 'min',   2,   'mean',  8.5,   'max',   16);
            map('d')   = struct('fixed', 0, 'min',   4,   'mean',  5.4,   'max',    8);
            map('fss') = struct('fixed', 0, 'min',   0.3, 'mean',  0.46,  'max',    0.6);
            map('p')   = struct('fixed', 0, 'min',   0.5, 'mean',  1.1,   'max',    1.5); 
            map('q0')  = struct('fixed', 0, 'min', 1e5,   'mean',  3.7e6, 'max',   50*max(this.dependentData));
            map('t0')  = struct('fixed', 0, 'min',   0,   'mean', 31,     'max', this.timeFinal/2); 

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('a'),  this.finalParams('d'), this.finalParams('fss'),this.finalParams('p'), ...
                this.finalParams('q0'), this.finalParams('t0'));
        end
        function Q    = estimateDataFast(this, a, d, fss, p, q0, t0)  
            Q = this.simulateQ(a, d, fss, p, q0, t0, this.timeInterpolants);
        end  
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end
