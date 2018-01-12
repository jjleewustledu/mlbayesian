classdef BretthorstMcmc < mlbayesian.IMcmcSolver & mlio.AbstractIO
	%% BRETTHORSTMCMC  

	%  $Revision$
 	%  was created 13-Dec-2017 20:00:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties        
        nBinsHist = 100 % nbins for hist
    end
    
	properties (Dependent)		
 		constructGenerative 
        isfinished
        kernel
        model
 	end

	methods
        
        %% GET/SET
        
        function g = get.constructGenerative(this)
            g = this.model.constructGenerative;
        end
        function g = get.isfinished(this)
            g = this.kernel.isfinished;
        end
        function g = get.kernel(this)
            g = this.kernel_;
        end
        function g = get.model(this)
            g = this.model_;
        end
        
        function this = set.constructGenerative(this, s)
            assert(islogical(s));
            assert(~isempty(this.model));
            this.model.constructGenerative = s;
        end
        function this = set.model(this, s)
            assert(isa(s, 'mlanalysis.IModel'));
            this.model_ = s;
            this.model_.classOfSolver = class(this);
            this.logger_ = mlpipeline.Logger(this.model_.fqfileprefix);
        end
        
        %% 
        
        function diagnose(this)
            this.kernel.fprintf;
            this.plotAnnealing;
            this.plotParameterCovariances;
            this.plotLogProbabilityQC;
            this.histStdOfError;
            this.model.plot;
            disp(this.model);
        end
        function this = estimateParameters(this)
            tic
            this.mcmcStruct_ = struct( ...
                'solverParameters', this.model_.solverParameters, ...
                'paramsPenalty', 0, ...
                'showAnnealing', true, ...
                'showBeta', true, ...
                'verbosity', 0, ...
                'adjustParams', @this.adjustParams__, ...
                'objectiveFunc', @this.objectiveFunc__, ...
                'parameterIndexToName', @this.parameterIndexToName__); 
            this.kernel_ = mlbayesian.McmcKernel.main(this.mcmcStruct_);
            this.model_ = this.model_.updateModel(this);
            toc
        end
        function plotAnnealing(this)
            %% PLOTANNEALING
            %  @returns subplots of \beta -> params
            %  @returns    plot  of \beta -> log(posterior)
            
            figure;
            krnl = this.kernel_;
            N = ceil(sqrt(double(krnl.nParams)));
            for p = 1:krnl.nParams                
                subplot(N, N, double(p));
                plot(1:krnl.nBeta, krnl.paramsBetas(p,:));
                xlabel('\beta (1/T)');
                ylabel(this.parameterIndexToName__(p));
            end
            
            figure;
            plot(1:krnl.nBeta, krnl.lpBetas);
            title('plotAnnealing');
            xlabel('\beta (1/T)');
            ylabel('log(posterior)');
        end
        function plotParameterCovariances(this)
            %% PCOLORPARAMETERCOVARIANCES
            %  @returns pcolor subplots of bivariate histograms of parameter covariances
            
            figure;
            krnl = this.kernel_;
            Np = double(krnl.nParams);
            for m = 1:Np % rows
                for n = 1:Np % cols
                    subplot(Np, Np, double(n + double(m-1)*Np));
                    hold on;
                    dat = [krnl.paramsHist(n,:)', krnl.paramsHist(m,:)'];
                    hn  = hist3(dat, [this.nBinsHist, this.nBinsHist]);
                    xb  = linspace(min(dat(:,1)), max(dat(:,1)), size(hn,1));
                    yb  = linspace(min(dat(:,2)), max(dat(:,2)), size(hn,1));
                    pcolor(xb, yb, hn);
                    shading('flat');
                    set(gca, 'YDir', 'Reverse');
                    xlabel(sprintf('%s x %s', ...
                        this.parameterIndexToName__(m), ...
                        this.parameterIndexToName__(n)));
                end
            end
        end
        function plotLogProbabilityQC(this)
            %% PLOTLOGPROBABILITYQC
            %  @returns overlaid plots of proposals -> log(prob(population))
            
            figure;
            hold on;
            krnl = this.kernel_;
            for k = 1:krnl.nPop
                plot(1:krnl.nProposalsQC, krnl.lpQC(k,:), 'Color', this.colorVariation(k, krnl.nPop));
            end
            hold off;
            title('plotLogProbabilityQC, populations in color');
            xlabel('proposal');
            ylabel('log(prob(population))');
        end
        function histStdOfError(this)
            figure;
            krnl = this.kernel_;
            hist(krnl.stdOfError, this.nBinsHist);
            title('histStdOfError');
            xlabel('std of error');
        end
        function save(this)
            save@mlio.AbstractIO;
            this.logger_.save;
        end
        function saveas(this, varargin)
            saveas@mlio.AbstractIO(varargin{:});
            this.logger_.saveas(varargin{:});
        end
		  
 		function this = BretthorstMcmc(varargin)
 			%% BRETTHORSTMCMC
            
            ip = inputParser;
            addParameter(ip, 'model', mlanalysis.NullModel, @(x) isa(x, 'mlanalysis.IModel'));
            addParameter(ip, 'datedFilename', true, @islogical);
            parse(ip, varargin{:});
            this.model = ip.Results.model;  
            this.logger_ = mlpipeline.Logger(this.model_.fqfileprefix, this);
            
            %% for mlio.AbstractIO
            this.fileprefix_ = strrep(class(this), '.', '_');
            if (ip.Results.datedFilename)
                this.fileprefix_ = [this.fileprefix_ '_' datestr(now, 30)];
            end
            this.filesuffix_ = '.mat';
 		end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        kernel_
        logger_
        mcmcStruct_
        model_
    end
    
    methods (Access = private)
        function varargout = adjustParams__(this, varargin)
            varargout{:} = this.model_.kernel.adjustParams(varargin{:});
        end
        function varargout = objectiveFunc__(this, varargin)
            varargout{:} = this.model_.kernel.objectiveFunc(varargin{:});
        end
        function varargout = parameterIndexToName__(this, varargin)
            varargout{:} = this.model_.parameterIndexToName(varargin{:});
        end
        function c = colorVariation(~, k, kmax)
            c   = [rvar(k/kmax) gvar(k/kmax) bvar(k/kmax)];
            
            function r = rvar(x)
                r = exp(-(x - 0.1667)^2/(0.0556)) + ...
                    exp(-(x - 1.1667)^2/(0.0556));
            end
            function g = gvar(x)
                g = exp(-(x - 0.5   )^2/(0.0556));
            end
            function b = bvar(x)
                b = exp(-(x - 0.8333)^2/(0.0556));
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
