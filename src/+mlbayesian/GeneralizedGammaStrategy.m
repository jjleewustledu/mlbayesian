classdef GeneralizedGammaStrategy < mlbayesian.AbstractMcmcStrategy
	%% GENERALIZEDGAMMASTRATEGY fits up to two generalized gamma distributions and a steady-state term to time-series data.
    %  https://en.wikipedia.org/wiki/Gamma_distribution#Characterization_using_shape_.CE.B1_and_rate_.CE.B2
    %  https://en.wikipedia.org/wiki/Generalized_gamma_distribution
    %  N.B.  \rho(t; a,b,p) = \frac{p b^a}{\Gamma(a/p)} t^{a-1} e^{-(b t)^p} \text{ with } 1/b > 0, a > 0, p > 0, t > 0.

    % Sample results:
    % Running mlbayesian_unittest.Test_GeneralizedGammaStrategy
    % annealing step 1 beta (1/temp) 2.040816e-02 logProb0 NaN
    % annealing step 2 beta (1/temp) 4.081633e-02 logProb0 3.629935e+00
    % annealing step 3 beta (1/temp) 6.122449e-02 logProb0 8.406637e+00
    % annealing step 4 beta (1/temp) 8.163265e-02 logProb0 3.098795e+01
    % annealing step 5 beta (1/temp) 1.020408e-01 logProb0 4.232267e+01
    % annealing step 6 beta (1/temp) 1.224490e-01 logProb0 5.567121e+01
    % annealing step 7 beta (1/temp) 1.428571e-01 logProb0 6.861425e+01
    % annealing step 8 beta (1/temp) 1.632653e-01 logProb0 7.801523e+01
    % annealing step 9 beta (1/temp) 1.836735e-01 logProb0 8.829775e+01
    % annealing step 10 beta (1/temp) 2.040816e-01 logProb0 1.045470e+02
    % 
    % annealing step 10 param  a1 mean 8.007840	 std 3.705102	 sigma 0.567130	 acc 0.114000
    % annealing step 10 param  b1 mean 0.340204	 std 0.109246	 sigma 0.023148	 acc 0.144000
    % annealing step 10 param  p1 mean 1.000000	 std 0.000000	 sigma 0.999000	 acc 0.000000
    % annealing step 10 param t01 mean 7.172316	 std 4.446294	 sigma 1.388889	 acc 0.163000
    % annealing step 10 param  a2 mean 19.819439	 std 8.663353	 sigma 2.835648	 acc 0.161000
    % annealing step 10 param  b2 mean 1.064188	 std 0.867900	 sigma 0.138889	 acc 0.205000
    % annealing step 10 param  p2 mean 1.000000	 std 0.000000	 sigma 0.999000	 acc 0.000000
    % annealing step 10 param t02 mean 2.827955	 std 2.441172	 sigma 4.018776	 acc 0.110000
    % annealing step 10 param weight1 mean 0.810204	 std 0.099160	 sigma 0.120563	 acc 0.145000
    % annealing step 10 param   S mean 0.000000	 std 0.000000	 sigma 0.100000	 acc 0.000000
    % annealing step 10 param   k mean 0.000000	 std 0.000000	 sigma 1.000000	 acc 0.000000
    % annealing step 10 param  t0 mean 40.000000	 std 0.000000	 sigma 10.000000	 acc 0.000000
    % annealing step 50 beta (1/temp) 1.020408e+00 logProb0 6.349310e+02
    % BEST-FIT    param  a1 value 4.187151
    % BEST-FIT    param  b1 value 0.181447
    % BEST-FIT    param  p1 value 1.000000
    % BEST-FIT    param t01 value 6.014813
    % BEST-FIT    param  a2 value 20.076335
    % BEST-FIT    param  b2 value 0.714215
    % BEST-FIT    param  p2 value 1.000000
    % BEST-FIT    param t02 value 1.044547
    % BEST-FIT    param weight1 value 0.582077
    % BEST-FIT    param   S value 0.000000
    % BEST-FIT    param   k value 0.000000
    % BEST-FIT    param  t0 value 40.000000
    
	%  $Revision$
 	%  was created 27-Jun-2017 22:11:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties 	
        a1  = 3.31
        b1  = 0.149
        p1  = 1
        t01 = 7.43 % t01 < t02 or t01 > t02 by this.adjustParams
        a2  = 20.2
        b2  = 0.700
        p2  = 1
        t02 = 0.130
        weight1 = 0.521 % > 0.5 by this.mapParams
        S   = 0
        k   = 0
        t0  = 40
        force_t01_lt_t02 = false
        
        notes = ''
        xLabel = 'time/s'
        yLabel = 'arbitrary'
    end
    
    properties (Dependent)
        baseTitle
        detailedTitle
        mapParams
    end

    methods (Static)
        function [this,lg] = doBayes
            this.filepath = fullfile(getenv('HOME'), 'Local', 'src', 'mlcvl', 'mlbayesian', 'data', '');
            cd(this.filepath);
            load('kernel6_span33_deg4.mat');
            %kernel = zeros(size(kernelBest));
            %kernel(12:46) = kernelBest(12:46);
            kernel = kernel(1:120);
            kernel = kernel/sum(kernel);
            this = mlbayesian.GeneralizedGammaStrategy({0:119}, {kernel'});
            [this,lg] = this.doItsBayes;
        end
        function        plotInitial
            this.filepath = fullfile(getenv('HOME'), 'Local', 'src', 'mlcvl', 'mlbayesian', 'data', '');
            cd(this.filepath);
            load('kernel6_span33_deg4.mat');
            %kernel = zeros(size(kernelBest));
            %kernel(12:46) = kernelBest(12:46);
            kernel = kernel(1:120); %#ok<NODEF>
            kernel = kernel/sum(kernel);
            this = mlbayesian.GeneralizedGammaStrategy({0:119}, {kernel'});
            this.plot;
        end
        function mdl  = model(varargin)
            mdl = mlbayesian.GeneralizedGammaStrategy.rho(varargin{:});
        end
        function r    = rho(a1, b1, p1, t01, a2, b2, p2, t02, weight1, S, k, t0, t)
            r = mlbayesian.GeneralizedGammaTerms.gammaStretchSeries( ...
                            a1, b1, p1, t01, a2, b2, p2, t02, weight1, t);
            %r = mlbayesian.GeneralizedGammaTerms.gammaPair(a1, b1, t01, t02, weight1, t);
        end
        function this = simulateMcmc(a1, b1, p1, t01, a2, b2, p2, t02, weight1, S, k, t0, t, mapParams, keysParams)
            import mlbayesian.*;
            rho = GeneralizedGammaStrategy.rho( ...
                a1, b1, p1, t01, a2, b2, p2, t02, weight1, S, k, t0, t);
            this = GeneralizedGammaStratey({t}, {rho});
            this.mapParams = mapParams;
            this.keysParams_ = keysParams;
            [this,lg] = this.doItsBayes;
            fprintf('%s\n', char(lg));
        end
    end
    
	methods 
        
        %% GET/SET
        
        function g    = get.baseTitle(this)
            g = sprintf('%s %s', class(this), pwd);
        end
        function dt   = get.detailedTitle(this)
            dt = sprintf('%s\na1 %g, b1 %g, p1 %g, t01 %g, \na2 %g, b2 %g, p2 %g, t02 %g, \nwt1 %g, S %g, k %g, t0 %g\n%S', ...
                         this.baseTitle, ...
                         this.a1, this.b1, this.p1, this.t01, ...
                         this.a2, this.b2, this.p2, this.t02, ...
                         this.weight1, this.S, this.k, this.t0, ...
                         this.notes);
        end
        function this = set.mapParams(this, m)
            assert(isa(m, 'containers.Map'));
            this.mapParams_ = m;
        end
        function m    = get.mapParams(this)
            if (~isempty(this.mapParams_))
                m = this.mapParams_;
                return
            end
            
            m = containers.Map;
            m('a1')  = struct('fixed', 0, 'min', 1,     'mean', this.a1,  'max', 50);
            m('b1')  = struct('fixed', 0, 'min', 0.000, 'mean', this.b1,  'max', 10);
            m('p1')  = struct('fixed', 1, 'min', 0.01,  'mean', this.p1,  'max', 10);
            m('t01') = struct('fixed', 0, 'min', 0,     'mean', this.t01, 'max', 100);
            m('a2')  = struct('fixed', 0, 'min', 1,     'mean', this.a2,  'max', 50);
            m('b2')  = struct('fixed', 0, 'min', 0.000, 'mean', this.b2,  'max', 10);
            m('p2')  = struct('fixed', 1, 'min', 0.01,  'mean', this.p2,  'max', 10);
            m('t02') = struct('fixed', 0, 'min', 0,     'mean', this.t02, 'max', 100);
            m('weight1') = struct('fixed', 0, 'min', 0.5+eps, 'mean', this.weight1, 'max', 1);
            m('S')   = struct('fixed', 1, 'min', 0,     'mean', this.S,   'max', 1);
            m('k')   = struct('fixed', 1, 'min', 0,     'mean', this.k,   'max', 10);
            m('t0')  = struct('fixed', 1, 'min', 0,     'mean', this.t0,  'max', 100);
        end
        
        %%
    
        function ps   = adjustParams(this, ps)
            theParams = this.theParameters;
            if (this.force_t01_lt_t02)                
                if (ps(theParams.paramsIndices('t01')) > ps(theParams.paramsIndices('t02')))
                    tmp                                = ps(theParams.paramsIndices('t01'));
                    ps(theParams.paramsIndices('t01')) = ps(theParams.paramsIndices('t02'));
                    ps(theParams.paramsIndices('t02')) = tmp;
                end
            else
                
                if (ps(theParams.paramsIndices('t02')) > ps(theParams.paramsIndices('t01')))
                    tmp                                = ps(theParams.paramsIndices('t02'));
                    ps(theParams.paramsIndices('t02')) = ps(theParams.paramsIndices('t01'));
                    ps(theParams.paramsIndices('t01')) = tmp;
                end
            end
        end  
        function ed   = estimateDataFast(this, a1, b1, p1, t01, a2, b2, p2, t02, weight1, S, k, t0)
            %% ESTIMATEDATAFAST is used by AbstractBayesianStrategy.theSolver.
            
            ed{1} = mlbayesian.GeneralizedGammaStrategy.rho(a1, b1, p1, t01, a2, b2, p2, t02, weight1, S, k, t0, this.times{1});
        end
        function this = estimateParameters(this, varargin)
            ip = inputParser;
            addOptional(ip, 'mapParams', this.mapParams, @(x) isa(x, 'containers.Map'));
            parse(ip, varargin{:});
            
            this = this.runMcmc(ip.Results.mapParams, 'keysToVerify', this.keysParams_);
        end
        function r    = itsRho(this)
            r =  mlbayesian.GeneralizedGammaStrategy.rho( ...
                this.a1, this.b1, this.p1, this.t01, this.a2, this.b2, this.p2, this.t02, this.weight1, this.S, this.k, this.t0, this.times{1});
        end
        function lg   = logging(this)
            lg = mlpipeline.Logger(this.fqfileprefix);
            lg.add('\n%s is working in %s\n', mfilename, pwd);
            if (~isempty(this.theSolver))
                s = this.theSolver;
                lg.add('bestFitParams -> %s\n', mat2str(s.bestFitParams));
                lg.add('meanParams -> %s\n', mat2str(s.meanParams));
                lg.add('stdParams -> %s\n', mat2str(s.stdParams));
            end
            lg.add('\n');
        end   
        function        plot(this, varargin)
            figure;
            plot(this.times{1}, this.itsRho, this.times{1}, this.dependentData{1}, varargin{:});
            title(this.detailedTitle, 'Interpreter', 'none');
            xlabel(this.xLabel);
            ylabel(this.yLabel);
        end  
        function        plotParVars(this, par, vars)
            assert(lstrfind(par, properties(this)));
            assert(isnumeric(vars));
            switch (par)
                case 'a1'
                    for v = 1:length(vars)
                        args{v} = { vars(v) this.b1 this.p1 this.t01 this.a2 this.b2 this.p2 this.t02 this.weight1 this.S  this.k  this.t0 };  %#ok<*AGROW>
                    end
                case 'b1'
                    for v = 1:length(vars)
                        args{v} = { this.a1 vars(v) this.p1 this.t01 this.a2 this.b2 this.p2 this.t02 this.weight1 this.S  this.k  this.t0 }; 
                    end
                case 'p1'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 vars(v) this.t01 this.a2 this.b2 this.p2 this.t02 this.weight1 this.S  this.k  this.t0 }; 
                    end
                case 't01'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 vars(v)  this.a2 this.b2 this.p2 this.t02 this.weight1 this.S  this.k  this.t0 }; 
                    end
                case 'a2'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 vars(v) this.b2 this.p2 this.t02 this.weight1 this.S  this.k  this.t0 };
                    end
                case 'b2'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 vars(v) this.p2 this.t02 this.weight1 this.S  this.k  this.t0 }; 
                    end
                case 'p2'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 this.b2 vars(v) this.t02 this.weight1 this.S  this.k  this.t0 }; 
                    end
                case 't02'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 this.b2 this.p2 vars(v)  this.weight1 this.S  this.k  this.t0 }; 
                    end
                case 'weight1'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 this.b2 this.p2 this.t02 vars(v)      this.S  this.k  this.t0 }; 
                    end
                case 'S'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 this.b2 this.p2 this.t02 this.weight1 vars(v) this.k  this.t0 }; 
                    end
                case 'k'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 this.b2 this.p2 this.t02 this.weight1 this.S  vars(v) this.t0 }; 
                    end
                case 't0'
                    for v = 1:length(vars)
                        args{v} = { this.a1 this.b1 this.p1 this.t01 this.a2 this.b2 this.p2 this.t02 this.weight1 this.S  this.k  vars(v) }; 
                    end
            end
            this.plotParArgs(par, args, vars);
        end
        function this = simulateItsMcmc(this)
            this = mlbayesian.GeneralizedGammaStratey.simulateMcmc( ...
                this.a1, this.b1, this.p1, this.t01, ...
                this.a2, this.b2, this.p2, this.t02, ...
                this.weight1, this.S, this.k, this.t0, this.times{1}, this.mapParams, this.keysParams_);
        end   
        function sse  = sumSquaredErrors(this, p)
            %% SUMSQUAREDERRORS returns the sum-of-square residuals for all cells of this.dependentData and 
            %  corresponding this.estimateDataFast.  Compared to AbstractMcmcStrategy.sumSquaredErrors, this 
            %  overriding implementation weights of the log-likelihood with Jeffrey's prior according to this.independentData.
            %  See also:  mlbayesian.AbstractMcmcStrategy.sumSquaredErrors, 
            %             mlkinetics.AbstractKinetics.jeffreysPrior.
            
            p   = num2cell(p);
            sse = 0;
            edf = this.estimateDataFast(p{:});
            for iidx = 1:length(this.dependentData)
                sse = sse + ...
                      sum( (this.dependentData{iidx} - edf{iidx}).^2 ); % .* ...
                            %mlbayesian.AbstractBayesianStrategy.slide( ...
                            %    this.jeffreysPrior{iidx}, this.independentData{iidx}, this.t01) );
            end
            if (sse < 10*eps)
                sse = sse + (1 + rand(1))*10*eps; 
            end
        end     
		  
 		function this = GeneralizedGammaStrategy(varargin)
 			%% GENERALIZEDGAMMASTRATEGY
 			%  @params independentData is cell.
            %  @params dependentData is cell.
            %  @returns this.

 			this = this@mlbayesian.AbstractMcmcStrategy(varargin{:}); 
            this.keysParams_ = {'a1' 'b1' 'p1' 't01' 'a2' 'b2' 'p2' 't02' 'weight1' 'S' 'k' 't0'};
 		end
 	end 
    
    %% PROTECTED
    
    properties (Access = 'protected')
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

