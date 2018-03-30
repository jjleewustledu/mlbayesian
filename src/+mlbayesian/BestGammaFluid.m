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
        baseTitle = 'BestGammaFluid'
        xLabel    = 'times/s'
        yLabel    = 'wellcounts'
        
        a  = 8.5
        d  = 5.4
        e  = 0.54
        p  = 1.1
        q0 = 3.7e6
        t0 = 31
    end 
    
    properties (Dependent)
        map
    end
    
    methods %% GET 
        function m = get.map(this)            
            m = containers.Map;
            m('a')  = struct('fixed', 0, 'min',   2,   'mean',  this.a,  'max',   16);
            m('d')  = struct('fixed', 0, 'min',   4,   'mean',  this.d,  'max',    8);
            m('e')  = struct('fixed', 0, 'min',   0.4, 'mean',  this.e,  'max',    0.7);
            m('p')  = struct('fixed', 0, 'min',   0.5, 'mean',  this.p,  'max',    1.5); 
            m('q0') = struct('fixed', 0, 'min', 1e5,   'mean',  this.q0, 'max',   50*max(this.dependentData));
            m('t0') = struct('fixed', 0, 'min',   0,   'mean',  this.t0, 'max', this.timeFinal/2); 
        end
    end
    
    methods (Static)
        function this = runGammaFluid(times, counts)
            this = mlbayesian.BestGammaFluid(times, counts);
            this = this.estimateParameters(this.map);
        end
        function Q    = fluidQ(a, d, e, p, q0, t0, times)
            idx_t0 = mlbayesian.BestGammaFluid.indexOf(times, t0);
            cnorm  = q0 * ((p/a^d)/gamma(d/p));
            exp1   = exp(-(times/a).^p);
            Q0     = cnorm * times.^(d-1) .* exp1;
            Q0     = abs(e * Q0 + (1 - e) * max(Q0) * (1 - exp1));
            
            Q             = zeros(1, length(times));
            Q(idx_t0:end) = Q0(1:end-idx_t0+1);
            assert(all(isreal(Q)), 'BestGammaFluid.fluidQ.Q was complex');
            assert(~any(isnan(Q)), 'BestGammaFluid.fluidQ.Q was NaN: %s', num2str(Q));
        end
        function this = simulateMcmc(a, d, fss, p, q0, t0, times, map)
            
            import mlbayesian.*;            
            Q    = BestGammaFluid.fluidQ(a, d, fss, p, q0, t0, times);
            this = BestGammaFluid(times, Q);
            this = this.estimateParameters(map) %#ok<NOPRT>
            
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
            this.expectedBestFitParams_ = [this.a this.d this.e this.p this.q0 this.t0]';
        end 
        function fq   = itsFluidQ(this)
            fq = mlbayesian.BestGammaFluid.fluidQ(this.a, this.d, this.e, this.p, this.q0, this.t0, this.times);
        end
        function this = estimateParameters(this, varargin)
            ip = inputParser;
            addOptional(ip, 'map', this.map, @(x) isa(x, 'containers.Map'));
            parse(ip, varargin{:});
            
            import mlbayesian.*;
            this.paramsManager = BayesianParameters(ip.Results.map);
            this.ensureKeyOrdering({'a' 'd' 'e' 'p' 'q0' 't0'});
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
            this.a  = this.finalParams('a');
            this.d  = this.finalParams('d');
            this.e  = this.finalParams('e');
            this.p  = this.finalParams('p');
            this.q0 = this.finalParams('q0');
            this.t0 = this.finalParams('t0');
        end
        function ed   = estimateData(this)
            keys = this.paramsManager.paramsMap.keys;
            ed = this.estimateDataFast( ...
                this.finalParams(keys{1}), ...
                this.finalParams(keys{2}), ...
                this.finalParams(keys{3}), ...
                this.finalParams(keys{4}), ...
                this.finalParams(keys{5}), ...
                this.finalParams(keys{6}));
        end
        function Q    = estimateDataFast(this, a, d, e, p, q0, t0)  
            Q = this.fluidQ(a, d, e, p, q0, t0, this.timeInterpolants);
        end  
        function x    = priorLow(~, x)
            x = 0.9 * x;
        end
        function x    = priorHigh(~, x)
            x = 1.1 * x;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end
