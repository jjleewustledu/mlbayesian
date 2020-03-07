classdef SumSquaredJeffreysResiduals < handle & mlbayesian.SumSquaredWeightedResiduals
	%% SUMSQUAREDJEFFREYSRESIDUALS  

	%  $Revision$
 	%  was created 28-Aug-2018 17:08:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    methods
        function plotJeffreysPrior(this)
            plot(this.jeffreysPrior_); 
            ylabel('probability');
            xlabel('time / s');            
            if (this.client_.uniformSampling)
                title('SumSquaredJeffreysResiduals.plotJeffreysPrior after buildJeffreysPrior');
            else
                title('SumSquaredJeffreysResiduals.plotJeffreysPrior after buildJeffreysPriorNonuniform');
            end
        end
        function e = sumSquaredErrors(this, p)
            %% SUMSQUAREDERRORS returns the sum-of-square residuals for all cells of this.client_.dependentData and 
            %  this.client_.estimateDataFast.  This implementation weights the residuals with Jeffrey's prior to
            %  describe a uniform prior probability for each decade of this.client_.independentData.
            p = num2cell(p);
            e = 0;
            edf = this.client_.estimateDataFast(p{:});
            for idd = 1:length(this.client_.dependentData)
                summand = abs(this.client_.dependentData{idd} - edf{idd}).^2 .* ...
                              this.client_.jeffreysPrior_{idd} ./ ...
                              this.client_.dependentData{idd};
                e = e + sum(summand(isfinite(summand)));
            end
            if (e < 10*eps)
                e = e + (1 + rand(1))*10*eps; 
            end
        end	
		  
 		function this = SumSquaredJeffreysResiduals(varargin)
 			%% SUMSQUAREDJEFFREYSRESIDUALS
 			%  @param .

 			this = this@mlbayesian.SumSquaredWeightedResiduals(varargin{:});
            
            if (this.client_.uniformSampling)
                this.jeffreysPrior_ = this.buildJeffreysPrior;
            else
                this.jeffreysPrior_ = this.buildJeffreysPriorNonuniform;                
            end
            this.plotJeffreysPrior;
 		end
    end
    
	methods (Access = protected)
        
        function tf = baseline(this)
            %  @return tf is cell-array.
            
            tf = cellfun(@(x) false(size(x)), this.client_.dependentData);
            for idd = 1:length(tf)
                signal     = this.client_.dependentData{idd};
                [~,imaxDD] = max( signal);
                meanDD     = mean(signal(1:imaxDD));
                stdDD      = std( signal(1:imaxDD));                
                tf{idd}(1:imaxDD) = signal(1:imaxDD) < min([meanDD abs(meanDD - stdDD)]);
            end
        end
        function p = buildJeffreysPrior(this)
            %% JEFFREYSPRIOR
            %  Cf. Gregory, Bayesian Logical Data Analysis for the Physical Sciences, sec. 3.7.1.
            %  @return p is cell-array.
            
            p = cell(this.client_.independentData);
            for iid = 1:length(p)
                t = this.client_.independentData{iid};
                for it = 1:length(t)
                    if (abs(t(it)) < eps)
                        t(it) = min(t(t > eps));
                    end
                end
                p{iid} = 1./(t*log(t(end)/t(1)));
                [~,endBaseline] = max(double(~this.baseline{iid}));
                p{iid}(1:endBaseline) = p{iid}(endBaseline);
                p{iid} = p{iid}/sum(p{iid});
                assert(all(p{iid} >= 0))
            end
        end
        function p = buildJeffreysPriorNonuniform(this)
            %% JEFFREYSPRIOR
            %  Cf. Gregory, Bayesian Logical Data Analysis for the Physical Sciences, sec. 3.7.1.
            %  @return p is cell-array.
            
            p = cell(this.client_.independentData);
            for iid = 1:length(p)
                t = this.client_.independentData{iid};
                for it = 1:length(t)
                    if (abs(t(it)) < eps)
                        t(it) = min(t(t > eps));
                    end
                end
                taus_ = t(2:end) - t(1:end-1);
                taus_ = [taus_ taus_(end)]; %#ok<AGROW>
                p{iid} = 1./(t*log(t(end)/t(1)));
                [~,endBaseline] = max(double(~this.baseline{iid}));
                p{iid}(1:endBaseline) = p{iid}(endBaseline);
                p{iid} = p{iid}.*taus_/sum(p{iid}.*taus_);
                assert(all(p{iid} >= 0))
            end
        end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        jeffreysPrior_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

