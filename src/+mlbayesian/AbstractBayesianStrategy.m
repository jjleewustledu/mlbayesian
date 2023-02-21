classdef (Abstract) AbstractBayesianStrategy < handle & mlio.AbstractHandleIO & mlbayesian.IBayesianStrategy
	%% ABSTRACTBAYESIANSTRATEGY  

	%  $Revision$
 	%  was created 23-Nov-2015 17:07:18
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  Copyright 2015-2017 John Joowon Lee.
 	
    
    properties (Dependent)
        
        % populate interp cells by setting base cells; interp uses NyquistFreqFactor
        
        baseTitle
        dependentData % cells, e.g., densities = f(time)
        dependentDataInterp % cells, interpolated by this.NyquistFreq
        dt % scalar := min(min(cell2mat(this.independentDeltas)))/this.NyquistFreqFactor
        independentData % cells, e.g., times
        independentDataInterp % cells, interpolated by this.NyquistFreq
        independentDeltas % repeats independentData(end) - independentData(end-1) so that 
                          % length(independentDeltas) = length(independentData).
        NyquistFreqFactor
        sessionData
        theParameters
        theSolver
        verbosity
        
        taus             % synonym of independentDeltas
        times            % synonym of independentData
        timeFinal        % independentData(end)
        timeInitial      % independentData(1)
        timeInterpolants % synonym of independentDataInterp
    end
    
    methods (Static)
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
        function h       = Heaviside(t, t0)
            h = zeros(size(t));
            h = h + double(t > t0);
        end
        function this    = load(fn) %#ok<STOUT>
            load(fn, 'this');
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
    
    methods 
        
        %% GET/SET
        
        function g = get.baseTitle(this)
            if (isempty(this.sessionData))
                g = sprintf('%s in %s', class(this), pwd);
                return
            end
            g = sprintf('%s in %s', class(this), this.sessionData.sessionFolder);
        end
        function g = get.dependentData(this)
            g = this.dependentData_;
        end
        function     set.dependentData(this, s)
            if (isempty(s))
                return
            end
            s = cellfun(@(x) this.ensureRow(x), s, 'UniformOutput', false);
            this.dependentData_ = s;
            for iis = 1:length(s)
                this.dependentDataInterp_{iis} = ...
                    pchip(this.independentData_{iis}, this.dependentData_{iis}, this.independentDataInterp{iis});                
            end
        end
        function g = get.dependentDataInterp(this)
            g  = this.dependentDataInterp_;
        end
        function g = get.dt(this)
            g = this.dt_;
        end
        function g = get.independentData(this)
            g = this.independentData_;
        end
        function     set.independentData(this, s)
            if (isempty(s))
                return
            end
            s = cellfun(@(x) this.ensureRow(x), s, 'UniformOutput', false);
            this.independentData_ = s;
            for iis = 1:length(s)
                this.independentDeltas_{iis} = ...
                    [(s{iis}(2:end) - s{iis}(1:end-1)) (s{iis}(end) - s{iis}(end-1))];
            end
            this.dt_ = min(min(cell2mat(this.independentDeltas_)))/this.NyquistFreqFactor;
            for iis = 1:length(s)
                this.independentDataInterp_{iis} = ...
                    this.independentData{iis}(1):this.dt_:this.independentData{iis}(end);                
            end
        end
        function g = get.independentDataInterp(this)
            g  = this.independentDataInterp_;
        end
        function g = get.independentDeltas(this)
            g = this.independentDeltas_;
        end
        function g = get.NyquistFreqFactor(this)
            g = this.NyquistFreqFactor_;
        end
        function g = get.sessionData(this)
            g = this.sessionData_;
        end
        function     set.sessionData(this, s)
            assert(isa(s, 'mlpipeline.ISessionData'));
            this.sessionData_ = s;
        end
        function g = get.theParameters(this)
            g = this.theParameters_;
        end
        function     set.theParameters(this, s)
            assert(isa(s, 'mlbayesian.IMcmcParameters'));
            this.theParameters_ = s;
        end
        function g = get.theSolver(this)
            g = this.theSolver_;
        end
        function     set.theSolver(this, s)
            assert(isa(s, 'mlbayesian.IMCMC'));
            this.theSolver_ = s;
        end
        function g = get.verbosity(this)
            g = this.verbosity_;
        end
        
        function g = get.taus(this)
            g = this.independentDeltas;
        end
        function g = get.times(this)
            g = this.independentData;
        end
        function g = get.timeFinal(this)
            for iidx = 1:length(this.independentData)
                g(iidx) = this.independentData{iidx}(end); %#ok<AGROW>
            end
        end
        function g = get.timeInitial(this) 
            for iidx = 1:length(this.independentData)
                g(iidx) = this.independentData{iidx}(1); %#ok<AGROW>
            end
        end
        function g = get.timeInterpolants(this)
            g = this.independentDataInterp;
        end
    
        %%
        
        function save(this)
            save(this.fqfilename, 'this');
        end
        
 		function this = AbstractBayesianStrategy(varargin)
 			%% ABSTRACTBAYESIANSTRATEGY
 			%  Usage:  this = AbstractBayesianStrategy([independent_data, dependent_data])
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addOptional(ip,  'indepData', {}, @iscell); 
            addOptional(ip,  'depData', {}, @iscell);
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.ISessionData'));
            addParameter(ip, 'mcmcParameters', [], @(x) isa(x, 'mlbayesian.IMcmcParameters'));
            addParameter(ip, 'NyquistFreqFactor', 2, @(x) isnumeric(x) && x > 1);
            parse(ip, varargin{:});            
 			
            this.NyquistFreqFactor_ = ip.Results.NyquistFreqFactor;
            this.independentData    = ip.Results.indepData;
            this.dependentData      = ip.Results.depData;
            this.sessionData_       = ip.Results.sessionData;
            this.theParameters_     = ip.Results.mcmcParameters;
            this.fileprefix         = sprintf(strrep(class(this), '.', '_'));
            this.filesuffix         = '.mat';
            this = this.setVerbosityCache;            
            for idd = 1:length(this.dependentData)
                assert(all(size(this.independentData{idd}) == size(this.dependentData{idd})));
            end
        end 
    end   
    
    %% PROTECTED
    
    properties (Access = 'protected')
        dependentData_
        dependentDataInterp_
        dt_
        independentData_
        independentDataInterp_
        independentDeltas_
        NyquistFreqFactor_
        sessionData_
        theParameters_
        theSolver_
        verbosity_
    end
    
    methods (Access = 'protected')
        function this = setVerbosityCache(this)
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

