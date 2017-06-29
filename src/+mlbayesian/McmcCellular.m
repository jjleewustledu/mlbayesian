classdef McmcCellular < mlbayesian.IMCMC 
	%% MCMCCELLULAR has the machinery to do a simple Bayesian estimation;
    %  it becomes more verbose for getenv('VERBOSITY') > 0 or for getenv('VERBOSE') == 1.

	%  $Revision$
 	%  was created 23-Nov-2015 17:37:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%  It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.  
    %% Copyright 2015 G. Larry Bretthorst, Joshua S. Shimony, John J. Lee.

    properties (Constant)
        NBINS       = 50   % nbins for hist
        FRACPEEK    =  0.2
        PARPEN      =  0.0 % -1.0 % minimal penalty for each param (unused)
        MAX_PROP    = 50
        FRAC_POPREP = 0.1
    end
    
    properties
        paramsBetas   
        paramsPopulations   
        paramsSigmas   
        paramsHist        
        annealingAvpar 
        annealingSdpar
        annealingInitz
        bestFitParams
        meanParams
        stdParams
        stdOfError 
        
        lpQC        
        lpBetas
        lpPopulations
        lpFinal        
    end
    
    properties (Dependent)
        parameters   
        nParams
        nProposals % = 100 number of loops in parameter prob phase
        nPop       % =  50 number of population
        nPopRep    %       number of population to replace
        nBeta      % =  50 number of temperature steps
        nAnneal    % =  20 number of loops per annealing temp
        nSamples        
        nProposalsQC
        showAnnealing
        showBeta
        showBestFit
        showFinalStats
        showPlots
        verbosity
    end
    
    methods %% GET/SET
        function n = get.parameters(this)
            n = this.mcmcStrategy_.theParameters;
        end
        function n = get.nParams(this)
            n = this.parameters.length;
        end
        function n = get.nProposals(this)
            n = this.parameters.nProposals;
        end
        function n = get.nPop(this)
            n = this.parameters.nPop;
        end
        function n = get.nPopRep(this)
            n = this.FRAC_POPREP * this.nPop;
        end
        function n = get.nBeta(this)
            n = this.parameters.nBeta;
        end
        function n = get.nAnneal(this)
            n = this.parameters.nAnneal;
        end
        function n = get.nSamples(this)
            n = this.parameters.nSamples;
        end
        function n = get.nProposalsQC(this)
            n = ceil(this.FRACPEEK * this.nProposals);
        end
        function tf = get.showAnnealing(this)
            tf = this.mcmcStrategy_.showAnnealing;
        end
        function tf = get.showBeta(this)
            tf = this.mcmcStrategy_.showBeta;
        end
        function tf = get.showBestFit(this)
            tf = this.mcmcStrategy_.showBestFit;
        end
        function tf = get.showFinalStats(this)
            tf = this.mcmcStrategy_.showFinalStats;
        end
        function tf = get.showPlots(this)
            tf = this.mcmcStrategy_.showPlots;
        end
        function v  = get.verbosity(this)
            v = this.mcmcStrategy_.verbosity;
        end
    end
    
	methods
        function this                = McmcCellular(mcmcStrat)
            
            p = inputParser;
            addRequired(p, 'mcmcStrat', @(x) isa(x, 'mlbayesian.IMcmcStrategy'));
            parse(p, mcmcStrat);
            
            this.mcmcStrategy_     = p.Results.mcmcStrat;    
            this.paramsBetas       = zeros(this.nParams, this.nBeta);  
            this.paramsPopulations = zeros(this.nParams, this.nPop);   
            this.paramsSigmas      = zeros(this.nParams, 1);  
            this.bestFitParams     = zeros(this.nParams, 1);  
            this.stdOfError        = zeros(this.nPop*this.nProposalsQC, 1); % follow the standard deviation of error     
            this.annealingAvpar    = zeros(this.nParams, 1);
            this.annealingSdpar    = zeros(this.nParams, 1);
            this.annealingInitz    = zeros(this.nParams, 1);
            
            this.lpBetas       = zeros(this.nBeta, 1);
            this.lpPopulations = zeros(this.nPop, 1);            
            this.lpQC          = zeros(this.nPop, this.nProposalsQC);              % qc on the lprob
            this.paramsHist    = zeros(this.nParams, this.nPop*this.nProposalsQC); % for histogram of parameters          
            
            %% %%%%%%%%%%%%%%%%%%%
            %% initialize the Mcmc
            %% %%%%%%%%%%%%%%%%%%%
            
            if (this.verbosity > eps)
                fprintf('mlbayesian.McmcCellular.ctor:  initializing McmcCellular'); end
            for m = 1:this.nPop
                this.paramsSigmas = (this.parameters.max - this.parameters.min)/10.0;
                for k = 1:this.nParams
                    if (this.parameters.fixed(k))
                        this.paramsPopulations(k,m) = this.parameters.fixedValue(k);
                    else
                        this.paramsPopulations(k,m) = this.parameters.mean(k) + this.parameters.std(k)*randn(1,1);
                        while (this.paramsPopulations(k,m)<this.parameters.min(k) || this.paramsPopulations(k,m)>this.parameters.max(k))
                            this.paramsPopulations(k,m) = this.parameters.mean(k) + this.parameters.std(k)*randn(1,1);
                        end
                    end
                end
            end
        end  
        function [parmax,avpar,this] = runMcmc(this)
            %% MCMCCELLULAR (Markov Chain Monte-Carlo) has the machinery to do a simple Bayesian estimation
            %
            %  after J. S. Shimony, Mar, 2014

            %% %%%%%%%%%%%%%%%%%%%%%%%
            %% annealing/burn-in phase
            %% %%%%%%%%%%%%%%%%%%%%%%%
            
            if (this.verbosity > eps)
                fprintf('mlbayesian.McmcCellular.runMcmc:  annealing/burn-in'); end
            lp0 = nan;
            for b = 1:this.nBeta  
                
                beta_ = (1.0/(this.nBeta-1.0))*b; 
                this.printBeta(b, beta_, lp0);
                parn = this.annealingInitz;
                this.lpBetas(b) = 0.0;

                %% population loop
                for m = 1:this.nPop
                    
                    ptmp = this.paramsPopulations(:,m);
                    lp0  = this.logProbability(ptmp, beta_, 1);
                    this.lpPopulations(m) = lp0;

                    for j = 1:this.nAnneal
                        for k = 1:this.nParams
                            
                            if (this.parameters.fixed(k)) 
                                continue; 
                            else
                                ptmp = this.paramsPopulations(:,m);
                                nprop = 0;
                                while (nprop < this.MAX_PROP)
                                    nprop = nprop + 1;
                                    dpar = this.paramsSigmas(k)*randn(1,1);
                                    if (ptmp(k)+dpar>=this.parameters.min(k) && ptmp(k)+dpar<=this.parameters.max(k)) 
                                        break; 
                                    end
                                end
                                if (nprop < this.MAX_PROP)
                                    ptmp(k) = ptmp(k) + dpar;
                                else
                                    if (1 == this.verbosity)
                                        fprintf('Warning: nprop %i, par %d, val %f, sigma %f\n',nprop,k,ptmp(k),this.paramsSigmas(k)); end
                                end
                                % new lp
                                [lp1,ptmp] = this.logProbability(ptmp, beta_, 1);

                                % metropolis acceptance
                                dele = lp1 - lp0;
                                if (dele > 0.0 || rand(1,1) < exp(dele))
                                    lp0 = lp1;
                                    this.lpPopulations(m) = lp1;
                                    parn(k) = parn(k) + 1;
                                    this.paramsPopulations(:,m) = ptmp;
                                end
                            end % end if (this.parameters.fixed(k))
                        end % end for k = 1:this.nParams
                    end % end for j = 1:this.nAnneal
                end % end for m = 1:this.nPop

                this = this.gatherAnnealingStats(b);
                this = this.BretthorstAdjustments(parn);    
                this.printAnnealing(b, parn);
                this = this.replacePoorMembers;

            end % end for b = 1:this.nBeta 

            if (this.showPlots); this.plotAnnealing; end
            
            %% %%%%%%%%%%%%%%%%%%%%%%%
            %% proposal/sampling phase
            %% %%%%%%%%%%%%%%%%%%%%%%%

            if (this.verbosity > eps)
                fprintf('mlbayesian.McmcCellular.runMcmc:  proposal/sampling'); end
            for m = 1:this.nPop
                
                ptmp = this.paramsPopulations(:,m);
                this.bestFitParams = this.paramsPopulations(:,m);
                [lp0,ptmp] = this.logProbability(ptmp, 1.0, 1);
                this.lpPopulations(m) = lp0;
                lpmax = lp0;
                parn = this.annealingInitz; 
                
                %% proposal (importance sampling) loop
                for j = 1:this.nProposals
                    for k = 1:this.nParams
                        
                        if (this.parameters.fixed(k)) 
                            continue; 
                        else
                            ptmp = this.paramsPopulations(:,m);
                            nprop = 0;
                            while (nprop < this.MAX_PROP)
                                nprop = nprop + 1;
                                dpar = this.paramsSigmas(k)*randn(1,1);
                                if (ptmp(k)+dpar>=this.parameters.min(k) && ptmp(k)+dpar<=this.parameters.max(k)) 
                                    break; 
                                end
                            end
                            if (nprop < this.MAX_PROP)
                                ptmp(k) = ptmp(k) + dpar;
                            else
                                if (1 == this.verbosity)
                                    fprintf('Warning: nprop %i, par %d, val %f, sigma %f\n',nprop,k,ptmp(k),this.paramsSigmas(k)); end
                            end
                            
                            % new lp
                            [lp1,ptmp] = this.logProbability(ptmp, 1.0, 1);

                            % store best for future use
                            if (lp1 > lpmax)
                                lpmax = lp1;
                                this.bestFitParams = ptmp;
                            end

                            % metropolis acceptance
                            dele = lp1 - lp0;
                            if (dele > 0.0 || rand(1,1) < exp(dele))
                                lp0 = lp1;
                                this.lpPopulations(m) = lp1;
                                parn(k) = parn(k) + 1;
                                this.paramsPopulations(:,m) = ptmp;
                            end
                        end % end if (this.parameters.fixed(k))
                    end % end for k = 1:this.nParams
                    
                    [lp1, this] = this.updateHistograms(j, m, lp0, lp1, ptmp);
                    
                end % end for j = 1:this.nProposals
            end % end for m = 1:this.nPop
            
            %% %%%%%%%%%
            %% reporting
            %% %%%%%%%%%
            
            parmax = this.bestFitParams;
            avpar  = this.annealingAvpar;
            this.lpFinal = lp1;

            this = this.printBestFit; 
            this = this.printFinalStats;
            if (this.showPlots)
                this.histParametersDistributions;
                this.plotLogProbabilityQC;
                this.histStdOfError;
            end
        end    
        function [lprob,paramsVec]   = logProbability(this, paramsVec, beta_, lpFlag)
            %% LOGPROBABILITY
            %
            %  lpFlag: -1  return Q
            %           0  return log of normalized Q
            %           1  return with beta and priors
            %
            %  from J. S. Shimony, Mar, 2014
            
            paramsVec = this.mcmcStrategy_.adjustParams(paramsVec);
            lprob = this.mcmcStrategy_.sumSquaredErrors(paramsVec);
            
            if (lpFlag == -1)
                return
            end

            % use log t-distribution for parameter degrees of freedom
            lprob = -0.5*this.nSamples*log(0.5*lprob);

            % add in beta and pretest probabilities
            if (lpFlag == 1)
                % beta only operates on likelihoods
                args  = this.PARPEN - (paramsVec - this.parameters.mean).^2 ./ (2 * this.parameters.std.^2);
                args  = double(~this.parameters.fixed) .* args;
                lprob = beta_ * lprob + sum(args);
            end
        end

        function this = printBestFit(this)
            % print best fit member of population
            
            if (~this.showBestFit); return; end
            for k = 1:this.nParams
                fprintf('BEST-FIT    param %3s value %f\n', this.paramIndexToLabel(k), this.bestFitParams(k));
            end
        end
        function this = printFinalStats(this)
            
            avpar = this.annealingInitz;
            sdpar = this.annealingInitz;
            for m = 1:this.nPop
                ptmp = this.paramsPopulations(:,m);
                for k = 1:this.nParams
                    avpar(k) = avpar(k) + ptmp(k);
                    sdpar(k) = sdpar(k) + ptmp(k)^2;
                end
            end
            for k = 1:this.nParams
                avpar(k) = avpar(k)/this.nPop;
                sdpar(k) = (sdpar(k) - this.nPop*avpar(k)^2)/(this.nPop - 1);
                if (sdpar(k) > 0.0)
                    sdpar(k) = sqrt(sdpar(k));
                else
                    sdpar(k) = 0.0;
                end
            end
            
            avpar = mean(this.paramsPopulations'); %#ok<UDIM>
            sdpar = std(this.paramsPopulations'); %#ok<UDIM>
            fprintf('\n');
            for k = 1:this.nParams
                this.meanParams(k) = avpar(k);
                this.stdParams(k)  = sdpar(k);
                if (this.showFinalStats)
                    fprintf('FINAL STATS param %3s mean  %f\t std %f\n', this.paramIndexToLabel(k), avpar(k), sdpar(k));
                end
            end
        end        
        
        function plotAnnealing(this)
            
            figure;            
            N = ceil(sqrt(double(this.nParams)));
            for k = 1:this.nParams                
                subplot(N, N, double(k));
                plot(1:this.nBeta, this.paramsBetas(k,:));
                xlabel(['plotAnnealing: ' this.paramIndexToLabel(k)]);
                ylabel('annealing paramsBetas');
            end
            
            figure;
            plot(1:this.nBeta, this.lpBetas);
            title('plotAnnealing');
            xlabel('beta (1/temp)');
            ylabel('log(prob)');
        end
        function histParametersDistributions(this)
            
            % histogram sampling phase
            % histogram parameter distribution
            figure;
            N = ceil(sqrt(double(this.nParams)));
            for k = 1:this.nParams
                subplot(N, N, double(k));
                hist(this.paramsHist(k,:), this.NBINS);
                xlabel(['Parameter ', this.paramIndexToLabel(k)]);
            end
        end
        function plotLogProbabilityQC(this)
            figure;
            hold on;
            for k = 1:this.nPop
                plot(1:this.nProposalsQC, this.lpQC(k,:), 'Color', this.colorVariation(k, this.nPop));
            end
            hold off;
            title('plotLogProbabilityQC');
            xlabel('proposal#');
            ylabel('log(prob)');
        end
        function histStdOfError(this)
            figure;
            hist(this.stdOfError, this.NBINS);
            title('histStdOfError');
            xlabel('std. dev. of error');
        end
    end
    
    %% PRIVATE
    
    properties (Access = 'private')
        mcmcStrategy_
    end
    
    methods (Access = 'private')
        function              printBeta(this, b, beta, lp0)
            if ( this.verbosity < 0.3); return; end
            if (~this.showBeta); return; end
            fprintf('annealing step %d beta (1/temp) %d logProb0 %d\n', b, beta, lp0); 
        end        
        function              printAnnealing(this, b, parn)
            if ( this.verbosity < 0.2); return; end
            if (~this.showAnnealing); return; end
            if ( ceil(this.FRACPEEK*this.nBeta) == b)
                fprintf('\n');
                for k = 1:this.nParams
                    fprintf('annealing step %d param %3s mean %f\t std %f\t sigma %f\t acc %f\n',...
                        b,this.paramIndexToLabel(k),this.annealingAvpar(k),this.annealingSdpar(k),this.paramsSigmas(k),parn(k)/(this.nAnneal*this.nPop));
                end
                fprintf('\n');
            end
        end
        
        function [this,lp1] = gatherAnnealingStats(this, b)
            this.annealingAvpar = this.annealingInitz;
            this.annealingSdpar = this.annealingInitz;
            for m = 1:this.nPop
                ptmp = this.paramsPopulations(:,m);
                for k = 1:this.nParams
                    this.annealingAvpar(k) = this.annealingAvpar(k) + ptmp(k);
                    this.annealingSdpar(k) = this.annealingSdpar(k) + ptmp(k)^2;
                end
                beta_ = (1.0/(this.nBeta-1.0))*b; 
                lp1 = this.logProbability(ptmp, beta_, 0);
                this.lpBetas(b) = this.lpBetas(b) + lp1/this.nPop;
            end
            for k = 1:this.nParams
                this.annealingAvpar(k) = this.annealingAvpar(k)/this.nPop;
                this.annealingSdpar(k) = (this.annealingSdpar(k) - this.nPop*this.annealingAvpar(k)^2)/(this.nPop - 1);
                if (this.annealingSdpar(k) > 0.0)
                    this.annealingSdpar(k) = sqrt(this.annealingSdpar(k));
                else
                    this.annealingSdpar(k) = 0.0;
                end
                this.paramsBetas(k,b) = this.annealingAvpar(k);
            end
        end
        function  this      = BretthorstAdjustments(this, parn)
            %% BRETTHORSTADJUSTMENTS adjust proposals al a Bretthorst

            for k = 1:this.nParams
                if (this.parameters.fixed(k)) 
                    continue;
                end
                
                if (parn(k) < 0.1*this.nAnneal*this.nPop)
                        this.paramsSigmas(k) = 0.2 *this.paramsSigmas(k);
                elseif (parn(k) < 0.2*this.nAnneal*this.nPop)
                    if (this.paramsSigmas(k) > 0.05*this.annealingSdpar(k))
                        this.paramsSigmas(k) =      this.paramsSigmas(k)/1.2;
                    end
                elseif (parn(k) > 0.8*this.nAnneal*this.nPop)
                    if (this.paramsSigmas(k) <     this.annealingSdpar(k))
                        this.paramsSigmas(k) = 5.0 *this.paramsSigmas(k);
                    end
                elseif (parn(k) > 0.6*this.nAnneal*this.nPop)
                    if (this.paramsSigmas(k) <     this.annealingSdpar(k))
                        this.paramsSigmas(k) = 2.0*this.paramsSigmas(k);
                    end
                elseif (parn(k) > 0.3*this.nAnneal*this.nPop)
                    if (this.paramsSigmas(k) <     this.annealingSdpar(k))
                        this.paramsSigmas(k) = 1.1*this.paramsSigmas(k);
                    end
                end
                if (this.paramsSigmas(k) > (this.parameters.max(k) - this.parameters.min(k)))
                    warning('mlbayesian:parameterOutOfBounds', 'McmcCellular.paramsSigmas(%i) too large\n', k);
                    this.paramsSigmas(k) =  this.parameters.max(k) - this.parameters.min(k);
                end           
            end  
        end
        function  this      = replacePoorMembers(this)
            for k = 1:this.nPopRep
                lpmax = -Inf;
                lpmin =  Inf;
                ilpmax = NaN;
                ilpmin = NaN;
                for m = 1:this.nPop
                    if (this.lpPopulations(m) == 0.0) 
                        continue; 
                    end
                    if (this.lpPopulations(m) < lpmin)
                        lpmin = this.lpPopulations(m);
                        ilpmin = m;
                    end
                    if (this.lpPopulations(m) > lpmax)
                        lpmax = this.lpPopulations(m);
                        ilpmax = m;
                    end
                end
                
                %% assure that these are not repeated
                if (isnan(ilpmin) || isnan(ilpmax))
                    continue
                end
                this.lpPopulations(ilpmin) = 0.0;
                this.lpPopulations(ilpmax) = 0.0;
                this.paramsPopulations(:,ilpmin) = this.paramsPopulations(:,ilpmax);
            end
        end
        function  c         = colorVariation(~, k, kmax)
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
        function  lbl       = paramIndexToLabel(this, idx)
            keys = this.parameters.keysParams;
            lbl  = keys{idx};
            if (isnumeric(lbl))
                lbl = num2str(lbl);
            end
        end
        function [lp1,this] = updateHistograms(this, j, m, lp0, lp1, ptmp)
            if (mod(j, ceil(1/this.FRACPEEK)) == 0)
                this.lpQC(m, j*this.FRACPEEK) = lp0;
                for k = 1:this.nParams
                    this.paramsHist(k,(m-1)*this.nProposalsQC + j*this.FRACPEEK) = ptmp(k);
                end
                lp1 = this.logProbability(ptmp, 1.0, -1);
                this.stdOfError((m-1)*this.nProposalsQC + j*this.FRACPEEK) = sqrt(lp1/(this.nSamples-2));
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

