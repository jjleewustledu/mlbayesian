classdef GeneralizedGammaProblem < mlbayesian.AbstractMcmcProblem
	%% GENERALIZEDGAMMAPROBLEM fits a generalized gamma to time-counts data.
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
        baseTitle = 'GeneralizedGammaProblem'
        xLabel    = 'times/s'
        yLabel    = 'wellcounts'
    end 
    
    methods (Static)  
        function this = simulateGG(a, d, p, q0, t0, times)
            %% SIMULATE a, d, p, t0
            
            import mlbayesian.*;
            Q    = GeneralizedGammaProblem.simulateQ(a, d, p, q0, t0, times);
            this = GeneralizedGammaProblem(times, Q);
            figure;
            plot(times, Q);
            title(sprintf('simulateGG:  a->%g, d->%g, p->%g, q0->%g, t0->%g, max(t)->%g', a, d, p, q0, t0, max(times)));
            xlabel('time/s');
            ylabel('generalized gamma/arbitrary');
        end
        function Q    = simulateQ(a, d, p, q0, t0, times)            
            idx_t0 = floor(t0) + 1;
            Q0     = q0 * (p/a^d) * times.^(d-1) .* exp(-(times/a).^p) / gamma(d/p); 
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
        end
        function this = simulateMcmc(a, d, p, q0, t0, times)
            this = mlbayesian.GeneralizedGammaProblem.simulateGG(a, d, p, q0, t0, times);
            this = this.estimateParameters;
        end
    end

	methods 		  
 		function this = GeneralizedGammaProblem(varargin) 
 			%% GENERALIZEDGAMMAPROBLEM 
 			%  Usage:  this = GeneralizedGammaProblem(times, well_counts) 

 			this = this@mlbayesian.AbstractMcmcProblem(varargin{:}); 
 		end 
        function this = estimateParameters(this)
            %% ESTIMATEPARAMETERS a, d, p, q0, t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('a')  = struct('fixed', 0, 'min', eps, 'mean', 10, 'max', 2560);
            map('d')  = struct('fixed', 0, 'min', eps, 'mean', 3,  'max', 768);
            map('p')  = struct('fixed', 0, 'min', eps, 'mean', 1,  'max', 256);  
            map('q0') = struct('fixed', 0, 'min', eps, 'mean', 1,  'max', 256); % max(this.dependentData)
            map('t0') = struct('fixed', 0, 'min',   0, 'mean', 10, 'max', 2560); % this.timeFinal  

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('a'), this.finalParams('d'), this.finalParams('p'), this.finalParams('q0'), this.finalParams('t0'));
        end
        function Q = estimateDataFast(this, a, d, p, q0, t0)  
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;            
            Q0     = q0 * (p/a^d) * times.^(d-1) .* exp(-(times/a).^p) / gamma(d/p);  
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
        end   
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

