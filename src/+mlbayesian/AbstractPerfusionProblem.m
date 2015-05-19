classdef (Abstract) AbstractPerfusionProblem < mlbayesian.AbstractMcmcProblem & mlbayesian.IPerfusionProblem 
	%% ABSTRACTPERFUSIONPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.5.0.197613 (R2015a) 
 	%  $Id$ 

    properties (Dependent)        
        concentration_a
        concentration_obs
        mtt_a
        mtt_obs
    end
    
    methods %% GET
        function ca = get.concentration_a(this)
            assert(~isempty(this.concentration_a_));
            ca = this.concentration_a_;
        end
        function cobs = get.concentration_obs(this)
            assert(~isempty(this.dependentData));
            cobs = this.dependentData;
        end
        function ca = get.mtt_a(this)
            assert(~isempty(this.mtt_a_));
            ca = this.mtt_a_;
        end
        function ca = get.mtt_obs(this)
            assert(~isempty(this.mtt_obs_));
            ca = this.mtt_obs_;
        end
    end
    
	methods (Static)        
        function this = load(varargin)  %#ok<VANUS>
            this = [];
        end        
        function ci   = concentration_i
            ci = [];
        end
        function args = interpolateData
            args = {};
        end
        
        function m = moment1(t, c)
            m = sum(t .* c) / sum(c);
        end  
        function [times,counts] = shiftData(times0, counts0, Dt)
            import mlbayesian.*
            if (Dt > 0)
                [times,counts] = AbstractPerfusionProblem.shiftDataRight(times0, counts0, Dt);
            else
                [times,counts] = AbstractPerfusionProblem.shiftDataLeft( times0, counts0, Dt);
            end
        end
        function [times,counts] = shiftDataLeft(times0, counts0, Dt)
            %  Dt in sec
            Dt     = abs(Dt);
            idx_0  = floor(sum(double(times0 < Dt + times0(1)))+1);
            times  = times0(idx_0:end);
            times  = times - times(1);
            counts = counts0(idx_0:end);
            counts = counts - min(counts);
        end
        function [times,counts] = shiftDataRight(times0, counts0, Dt)
            %  Dt in sec
            Dt     = abs(Dt);
            lenDt  = ceil(Dt/(times0(2) - times0(1)));
            newLen = length(counts0) + lenDt;
            
            times0 = times0 - times0(1) + Dt;
            times  = [0:1:lenDt-1 times0];
            counts = counts0(1) * ones(1,newLen);            
            counts(end-length(counts0)+1:end) = counts0;
            counts = counts - min(counts);
        end
        function this = simulateMcmc
            this = [];
        end
    end 
    
    methods
        function this = AbstractPerfusionProblem(conc_a, varargin)
 			this = this@mlbayesian.AbstractMcmcProblem(varargin{:});
            assert(isnumeric(conc_a));
            this.concentration_a_ = conc_a;
        end
        function this = estimateAll(this)
            this = this.estimateParameters(this.map);
        end        
        function this = estimateParameters(this)
        end
        function ed   = estimateData(this) %#ok<MANU>
            ed = [];
        end
        function ed   = estimateDataFast(this) %#ok<MANU>
            ed = [];
        end
        function ci = itsConcentration_i(this) %#ok<MANU>
            ci = [];
        end
        function this = save(this)
            this = this.saveas(sprintf('%s.save.mat', class(this)));
        end
        function this = saveas(this, fn)
            abstractPerfusionProblem = this; %#ok<NASGU>
            save(fn, 'abstractPerfusionProblem');   
        end
        function this = simulateItsMcmc(this) %#ok<MANU>
            this = [];
        end        
        
        function plotInitialData(this)
            figure;
            max_a   = max(this.concentration_a_);
            max_obs = max(this.concentration_obs);
            plot(this.times, this.concentration_a_/max_a, ...
                 this.times, this.concentration_obs/max_obs);
            title(sprintf('%s plotInitialData', this.baseTitle), 'Interpreter', 'none');
            legend('concentration_a', 'concentration_{obs}');
            xlabel(this.xLabel);
            ylabel(sprintf('%s; rescaled %g, %g', this.yLabel, max_a, max_obs));
        end        
        function plotProduct(this)
            figure;
            plot(this.times, this.itsConcentration_i, this.times, this.concentration_obs, 'o');
            legend('concentration_i', 'concentration_{obs}');
            title(this.detailedTitle, 'Interpreter', 'none');
            xlabel(this.xLabel);
            ylabel(this.yLabel);
        end
        function plotParVars(this) %#ok<MANU>
        end
    end
    
    %% PROTECTED
    
    properties (Access = 'protected')
        concentration_a_
        mtt_a_
        mtt_obs_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

