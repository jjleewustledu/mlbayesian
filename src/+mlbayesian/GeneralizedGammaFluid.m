classdef GeneralizedGammaFluid < mlbayesian.AbstractMcmcProblem
	%% GENERALIZEDGAMMAFLUID fits a generalized gamma & exponential recovery to time-counts data.
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
        baseTitle = 'GeneralizedGammaFluid'
        xLabel    = 'times/s'
        yLabel    = 'wellcounts'
    end 
    
    methods (Static)  
        function Q    = simulateQ(a, d, da, dp, fss, p, q0, t0, times)            
            idx_t0 = floor(t0) + 1;            
            cnorm  = q0 * ((p/a^d)/gamma(d/p));
            Q0     = cnorm * times.^(d-1) .* exp(-(times/a).^p);
            Q0     = Q0 + ...
                        max(Q0) * fss * (1 - exp(-(times/(a+da)).^(p+dp))); 
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
        end
        function this = simulateMcmc(a, d, da, dp, fss, p, q0, t0, times)
            
            import mlbayesian.*;            
            Q    = GeneralizedGammaFluid.simulateQ(a, d, da, dp, fss, p, q0, t0, times);
            this = GeneralizedGammaFluid(times, Q);
            this = this.estimateParameters %#ok<NOPRT>
            
            figure;
            plot(times, Q, ...
                 times, GeneralizedGammaFluid.simulateQ( ...
                        this.finalParams('a'), ...
                        this.finalParams('d'), ...
                        this.finalParams('da'), ...
                        this.finalParams('dp'), ...
                        this.finalParams('fss'), ...
                        this.finalParams('p'), ...
                        this.finalParams('q0'), ...
                        this.finalParams('t0'), times), 'o');
            title(sprintf('simulateMcmc expected:  a->%g, d->%g, da->%g, dp->%g, fss->%g, p->%g, q0->%g, t0->%g, max(t)->%g', ...
                  a, d, da, dp, fss, p, q0, t0, max(times)));
            xlabel('time/s');
            ylabel('generalized gamma/arbitrary');
        end
    end

	methods 		  
 		function this = GeneralizedGammaFluid(varargin) 
 			%% GeneralizedGammaFluid 
 			%  Usage:  this = GeneralizedGammaFluid(times, well_counts) 

 			this = this@mlbayesian.AbstractMcmcProblem(varargin{:}); 
 		end 
        function this = estimateParameters(this)
            %% ESTIMATEPARAMETERS a, d, p, q0, t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('a')   = struct('fixed', 0, 'min', eps, 'mean', 10,   'max', this.timeFinal);
            map('d')   = struct('fixed', 0, 'min', eps, 'mean',  3,   'max',    8);
            map('fss') = struct('fixed', 0, 'min',   0, 'mean',  0.1, 'max',    0.5);
            map('p')   = struct('fixed', 0, 'min', eps, 'mean',  1,   'max',    8); 
            map('q0')  = struct('fixed', 0, 'min', eps, 'mean',  1,   'max',  1e2*max(this.dependentData));
            map('t0')  = struct('fixed', 0, 'min',   0, 'mean', 10,   'max', this.timeFinal);            
            
            maxa = map('a').max;
            maxp = map('p').max;
            meanfss = map('fss').mean;
            map('da')  = struct('fixed', 0, 'min', -maxa*meanfss, 'mean', 0, 'max', maxa*meanfss);
            map('dp')  = struct('fixed', 0, 'min', -maxp*meanfss, 'mean', 0, 'max', maxp*meanfss);

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('a'),  this.finalParams('d'), this.finalParams('da'), this.finalParams('dp'),  ...
                this.finalParams('fss'),this.finalParams('p'), this.finalParams('q0'), this.finalParams('t0'));
        end
        function Q    = estimateDataFast(this, a, d, da, dp, fss, p, q0, t0)  
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            cnorm  = q0 * ((p/a^d)/gamma(d/p));
            Q0     = cnorm * times.^(d-1) .* exp(-(times/a).^p);
            Q0     = Q0 + ...
                        max(Q0) * fss * (1 - exp(-(times/(a+da)).^(p+dp)));    
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
        end   
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

