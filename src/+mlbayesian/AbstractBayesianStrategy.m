classdef (Abstract) AbstractBayesianStrategy < mlbayesian.IBayesianStrategy
	%% ABSTRACTBAYESIANSTRATEGY  

	%  $Revision$
 	%  was created 23-Nov-2015 17:07:18
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	
    properties         
        independentData % cells, e.g., times
        dependentData   % cells, e.g., densities = f(time)
        theParameters
        theSolver
    end
    
    properties (Dependent)
        expectedBestFitParams
        bestFitParams
        meanParams
        stdParams
        stdOfError
        verbosity
    end
    
    methods %% GET
        function e = get.expectedBestFitParams(this)
            assert(~isempty(this.expectedBestFitParams_), ...
                   'mlbayesian:attemptToAccessUnassignedVar', ...
                   'concrete implementation of AbstractBayesianStrategy must assign this.expectedBestFitParams_');
            e = this.expectedBestFitParams_;
        end
        function p = get.bestFitParams(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.bestFitParams;
        end
        function p = get.meanParams(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.meanParams;
        end
        function p = get.stdParams(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.stdParams;
        end
        function p = get.stdOfError(this)
            assert(~isempty(this.theSolver));
            p = this.theSolver.stdOfError;
        end
        function v = get.verbosity(this)
            v = this.verbosity_;
        end
    end
    
	methods 
 		function this = AbstractBayesianStrategy(varargin)
 			%% ABSTRACTBAYESIANSTRATEGY
 			%  Usage:  this = AbstractBayesianStrategy([independent_data, dependent_data])
            
            p = inputParser;
            addOptional(p, 'indepData', [], @iscell); 
            addOptional(p,   'depData', [], @iscell);
            parse(p, varargin{:});            
 			
            this.independentData = p.Results.indepData;
            this.dependentData   = p.Results.depData;
            for didx = 1:length(this.dependentData)
                assert(all(size(this.independentData{didx}) == size(this.dependentData{didx})));
            end
            this = this.setVerbosityCache;
        end 
    end
    
    methods (Static)
        function h    = Heaviside(t, t0)
            h = zeros(size(t));
            h = h + double(t > t0);
        end
        function [vec,T] = ensureRow(vec)
            if (~isrow(vec))
                vec = vec';
                T = true;
                return
            end
            T = false;
        end
        function tf   = uniformSampling(t)
            t   = mlsystem.VectorTools.ensureRowVector(t);
            dts = t(2:end) - t(1:end-1);
            dt1 = t(2) - t(1);
            tf  = all(abs(dt1*ones(1,length(dts)) - dts) < eps('single'));
        end
        function conc = slide(conc, t, Dt)
            %% SLIDE works for inhomogeneous t according to the ability of pchip to interpolate
            
            import mlbayesian.AbstractBayesianStrategy.*;
            [conc,trans] = ensureRow(conc);
            t          = ensureRow(t);
            
            tspan = t(end) - t(1);
            tinc  = t(2) - t(1);
            t_    = [(t - tspan - tinc) t];   % prepend times
            conc_ = [zeros(size(conc)) conc]; % prepend zeros
            conc  = pchip(t_, conc_, t - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts to right
            
            if (trans)
                conc = conc';
            end
        end
    end    
    
    %% PROTECTED
    
    properties (Access = 'protected')
        expectedBestFitParams_
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

