classdef AbstractDynamicProblem < mlbayesian.AbstractMcmcProblem
	%% ABSTRACTDYNAMICPROBLEM contains empty methods to establish abstractions. 

	%  $Revision$
 	%  was created 25-Jan-2017 18:17:56
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlsiemens/src/+mlsiemens.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    properties
        dependentData1
        dependentData11
    end
    
    properties (Dependent)
        manifold
    end
    
    methods % GET
        function g = get.manifold(this)
            g = this.manifold_;
        end
    end
    
	methods (Static)    
        function d              = estimatedData
            d = [];
        end
        function args           = interpolateData
            args = {};
        end
        function this           = load(varargin)  %#ok<VANUS>
            this = [];
        end            
        function tito           = indexTakeOff(curve)
            maxCurve = max(curve);
            minCurve = min(curve);
            for ti = 1:length(curve)
                if (curve(ti) - minCurve > 0.05 * (maxCurve - minCurve))
                    break;
                end
            end
            tito = max(1, ti - 1);
        end
        function m              = moment1(t, c)
            import mlbayesian.*;
            tto = t(AbstractDynamicProblem.indexTakeOff(c));
            m = sum((t - tto) .* c) / sum(c);
        end
        function c              = myPchip(t0, c0, t)
            if (t(end) > t0(end))
                t0(end+1) = t(end);
                c0(end+1) = c0(end);
            end
            c = pchip(t0, c0, t);
        end
        function [times,counts] = shiftData(times0, counts0, Dt)
            import mlbayesian.*
            if (Dt > 0)
                [times,counts] = AbstractDynamicProblem.shiftDataRight(times0, counts0, Dt);
            else
                [times,counts] = AbstractDynamicProblem.shiftDataLeft( times0, counts0, Dt);
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
        function this           = simulateMcmc
            this = [];
        end
    end 
    
    methods
        function this = AbstractDynamicProblem(varargin)
            ip = inputParser;
            addRequired(ip, 'times', @isnumeric);
            addRequired(ip, 'manifold', @iscell);
            parse(ip, varargin{:});
            
            this.manifold_ = ip.Results.manifold;
            this.independentData = ip.Results.times;
            this.dependentData   = this.manifold{1};
            this.dependentData1  = this.manifold{2};
            this.dependentData11 = this.manifold{3};
        end
        function this = estimateAll(this)
            this = this.estimateParameters(this.map);
        end        
        function this = estimateParameters(this)
        end
        function ed   = estimateData(~)
            ed = [];
        end
        function ed   = estimateDataFast(~)
            ed = [];
        end
        function d    = itsEstimatedData(~)
            d = [];
        end 
        function d1   = itsEstimatedData1(~)
            d1 = [];
        end
        function d11  = itsEstimatedData11(~)
            d11 = [];
        end
        function this = save(this)
            this = this.saveas(sprintf('%s.save.mat', class(this)));
        end
        function this = saveas(this, fn)
            abstractDynamicProblem = this; %#ok<NASGU>
            save(fn, 'abstractDynamicProblem'); % matlab-system save the object abstractDynamicProblem
        end
        function this = simulateItsMcmc(this) %#ok<MANU>
            this = [];
        end        
        
        function plotInitialData(this)
            figure;
            max_y1 = max(this.dependentData1);
            max_y  = max(this.dependentData);
            plot(this.times, this.dependentData1/max_y1, ...
                 this.times, this.dependentData/max_y);
            title(sprintf('%s plotInitialData', this.baseTitle), 'Interpreter', 'none');
            legend('dependentData1', 'dependentData');
            xlabel(this.xLabel);
            ylabel(sprintf('%s; rescaled %g, %g', this.yLabel, max_y1, max_y));
        end        
        function plotProduct(this)
            figure;
            plot(this.times, this.dependentData1, this.times, this.dependentData, 'o');
            legend('dependentData1', 'dependentData');
            title(this.detailedTitle, 'Interpreter', 'none');
            xlabel(this.xLabel);
            ylabel(this.yLabel);
        end
        function plotParVars(this) %#ok<MANU>
        end
    end
    
    %% PROTECTED
    
    properties (Access = 'protected')
        manifold_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

