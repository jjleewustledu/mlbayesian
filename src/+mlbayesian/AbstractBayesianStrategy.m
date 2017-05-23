classdef (Abstract) AbstractBayesianStrategy < mlio.AbstractIO & mlbayesian.IBayesianStrategy
	%% ABSTRACTBAYESIANSTRATEGY  

	%  $Revision$
 	%  was created 23-Nov-2015 17:07:18
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  Copyright 2015-2017 John Joowon Lee.
 	
    
    properties (Dependent)
        bestFitParams
        dependentData   % cells, e.g., densities = f(time)
        expectedBestFitParams
        independentData % cells, e.g., times
        meanParams
        sessionData
        stdParams
        stdOfError
        theParameters
        theSolver
        verbosity
    end
    
    methods %% GET/SET
        function p    = get.bestFitParams(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.bestFitParams;
        end
        function g    = get.dependentData(this)
            assert(~isempty(this.dependentData_));
            g = this.dependentData_;
        end
        function this = set.dependentData(this, s)
            assert(iscell(s));
            this.dependentData_ = s;
        end
        function e    = get.expectedBestFitParams(this)
            assert(~isempty(this.expectedBestFitParams_), ...
                   'mlbayesian:attemptToAccessUnassignedVar', ...
                   'concrete implementation of AbstractBayesianStrategy must assign this.expectedBestFitParams_');
            e = this.expectedBestFitParams_;
        end
        function g    = get.independentData(this)
            assert(~isempty(this.independentData_));
            g = this.independentData_;
        end
        function this = set.independentData(this, s)
            assert(iscell(s));
            this.independentData_ = s;
        end
        function p    = get.meanParams(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.meanParams;
        end
        function g    = get.sessionData(this)
            g = this.sessionData_;
        end
        function this = set.sessionData(this, s)
            assert(isa(s, 'mlpipeline.ISessionData'));
            this.sessionData_ = s;
        end
        function p    = get.stdParams(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.stdParams;
        end
        function p    = get.stdOfError(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.stdOfError;
        end
        function g    = get.theParameters(this)
            assert(~isempty(this.theParameters_));
            g = this.theParameters_;
        end
        function this = set.theParameters(this, s)
            assert(isa(s, 'mlbayesian.IMcmcParameters'));
            this.theParameters_ = s;
        end
        function g    = get.theSolver(this)
            g = this.theSolver_;
        end
        function this = set.theSolver(this, s)
            assert(isa(s, 'mlbayesian.IMCMC'));
            this.theSolver_ = s;
        end
        function v    = get.verbosity(this)
            v = this.verbosity_;
        end
    end
    
    methods (Static)
        function this = load(fn) %#ok<STOUT>
            load(fn, 'this');
        end
    end
	methods 
 		function this = AbstractBayesianStrategy(varargin)
 			%% ABSTRACTBAYESIANSTRATEGY
 			%  Usage:  this = AbstractBayesianStrategy([independent_data, dependent_data])
            
            p = inputParser;
            addOptional(p, 'indepData', {[]}, @iscell); 
            addOptional(p,   'depData', {[]}, @iscell);
            parse(p, varargin{:});            
 			
            this.independentData = p.Results.indepData;
            this.dependentData   = p.Results.depData;
            for didx = 1:length(this.dependentData)
                assert(all(size(this.independentData{didx}) == size(this.dependentData{didx})));
            end
            this = this.setVerbosityCache;            
            this.fileprefix = sprintf(strrep(class(this), '.', '_'));
            this.filesuffix = '.mat';
        end 
        function save(this)
            save(this.fqfilename, 'this');
        end
    end
    
    methods (Static)
        function h       = Heaviside(t, t0)
            h = zeros(size(t));
            h = h + double(t > t0);
        end
        function           appendState(fname, state)
            if (~lstrfind(fname, '.mat'))
                fname = [fname '.mat'];
            end
            state1 = state;
            fields1 = fieldnames(state1);
            load(fname, 'state');
            for f = 1:length(fields1)
                if (~isfield(state, fields1{f}))
                    state.(fields1{f}) = state1.(fields1{f});
                end
            end
            save(fname, 'state');
        end
        function [vec,T] = ensureRow(vec)
            if (~isrow(vec))
                vec = vec';
                T = true;
                return
            end
            T = false; 
        end
        function [t,interp1,interp2] = interpolateAll(t1, conc1, t2, conc2)
            %% INTERPOLATEALL interpolates variably sampled {t1 conc1} and {t2 conc2} to {t interp1} and {t interp2}
            %  so that t satisfies Nyquist sampling
            
            dt = min([timeDifferences(t1) timeDifferences(t2)]) / 2;
            tInf = min([t1 t2]);
            tSup = max([t1 t2]);
            
            t  = tInf:dt:tSup;
            interp1 = pchip(t1,conc1,t);
            interp2 = pchip(t2,conc2,t);

            function timeDiffs = timeDifferences(times)
                timeDiffs = times(2:end) - times(1:end-1);
            end
        end
        function           saveState(fname, state)
            if (isempty(fname) || isempty(state))
                return
            end
            if (~lstrfind(fname, '.mat'))
                fname = [fname '.mat'];
            end
            save(fname, 'state');
        end
        function conc    = slide(conc, t, Dt)
            %% SLIDE slides discretized function conc(t) to conc(t - Dt);
            %  Dt > 0 will slide conc(t) towards later times t.
            %  Dt < 0 will slide conc(t) towards earlier times t.
            %  It works for inhomogeneous t according to the ability of pchip to interpolate.
            %  It may not preserve information according to the Nyquist-Shannon theorem.  
            
            import mlbayesian.*;
            [conc,trans] = AbstractBayesianStrategy.ensureRow(conc);
            t            = AbstractBayesianStrategy.ensureRow(t);
            
            tspan = t(end) - t(1);
            tinc  = t(2) - t(1);
            t_    = [(t - tspan - tinc) t];   % prepend times
            conc_ = [zeros(size(conc)) conc]; % prepend zeros
            conc_(isnan(conc_)) = 0;
            conc  = pchip(t_, conc_, t - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts to right
            
            if (trans)
                conc = conc';
            end
        end
        function A       = pchip(t, A, t_, Dt)
            %% PCHIP slides discretized function A(t) to A(t_ - Dt);
            %  Dt > 0 will slide conc(t) towards to later values of t.
            %  Dt < 0 will slide conc(t) towards to earlier values of t.
            %  It works for inhomogeneous t according to the ability of pchip to interpolate.
            %  It may not preserve information according to the Nyquist-Shannon theorem.  
            %  @param t  is the initial t sampling
            %  @param A  is the initial A sampling
            %  @param t_ is the final   t sampling
            %  @param Dt is the shift of t_
            
            tspan = t(end) - t(1);
            dt    = t(2) - t(1);
            t     = [(t - tspan - dt) t]; % prepend times
            A     = [zeros(size(A)) A]; % prepend zeros
            A     = pchip(t, A, t_ - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts conc to right
        end
        function tf      = uniformSampling(t)
            t   = mlsystem.VectorTools.ensureRowVector(t);
            dts = t(2:end) - t(1:end-1);
            dt1 = t(2) - t(1);
            tf  = all(abs(dt1*ones(1,length(dts)) - dts) < eps('single'));
        end
    end    
    
    %% PROTECTED
    
    properties (Access = 'protected')
        dependentData_
        expectedBestFitParams_
        independentData_
        sessionData_
        theParameters_
        theSolver_
        verbosity_
    end
    
    methods (Access = 'protected')
        function this =       setVerbosityCache(this)
            this.verbosity_ = str2num(getenv('VERBOSITY')); %#ok<ST2NM>
            if (isempty(this.verbosity_))
                this.verbosity_ = str2num(getenv('VERBOSE'));  %#ok<ST2NM>
            end 
            if (~isempty(this.verbosity_))
                fprintf('mlbayesian.AbstractBayesianStrategy.setVerbosityCache:\n');
                fprintf('\tverbosity is %g;\n', this.verbosity_); 
                fprintf('\tadjust by setting ENV variable VERBOSITY to [0 1] or VERBOSE to true/false.\n');
            end
        end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

