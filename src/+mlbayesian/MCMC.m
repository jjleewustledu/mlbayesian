classdef MCMC < mlbayesian.IMCMC 
	%% MCMC has the machinery to do a simple Bayesian estimation 

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 

    properties (Constant)
        NPROPOSALS = 100   % number of loops in parameter prob phase
        NPOP       =  50   % number of population
        NPOPREP    =   5   % number of population to replace
        NBETA      =  50   % number of temperature steps
        NANNEAL    =  20   % number of loops per annealing temp
        NMOD       =   1   % number of models (unused)
        NBINS      =  50   % nbins for hist
        FRACPEEK   =   0.2
        PARPEN     =   0.0 % -1.0 % minimal penalty for each param (unused)
        LRG        =   1.0e20
    end
    
    properties
        dependentData        
        paramsData   
        paramsBetas   
        paramsPopulations   
        paramsSigmas          
        annealingAvpar 
        annealingSdpar
        annealingInitz
        bestFitParams
        
        lpBetas
        lpPopulations
        lpFinal
        
        paramsHist 
        logProbQC  
        stdOfError 
        
        showBestFit    = true
        showFinalStats = true      
        verbosity      = 0
    end
    
    properties (Dependent)
        nParams
        nSamples
        nProposalsQC
        showPlots
    end
    
    methods %% GET/SET
        function n = get.nParams(this)
            n = this.paramsData.length;
        end
        function n = get.nSamples(this)
            n = length(this.dependentData);
        end
        function n = get.nProposalsQC(this)
            n = this.FRACPEEK * this.NPROPOSALS;
        end
        function tf = get.showPlots(this)
            tf = this.bayesianProblem_.showPlots;
        end
    end
    
	methods        
        function this                = MCMC(bayesProb, depDat, paramsDat)
            
            p = inputParser;
            addRequired(p, 'bayesProb', @(x) isa(x, 'mlbayesian.IBayesianProblem'));
            addRequired(p, 'depDat',    @isnumeric);
            addRequired(p, 'paramsDat', @(x) isa(x, 'mlbayesian.IBayesianParameters'));
            parse(p, bayesProb, depDat, paramsDat);            
            
            this.bayesianProblem_  = p.Results.bayesProb;
            this.dependentData     = p.Results.depDat;
            this.paramsData        = p.Results.paramsDat;   
            this.paramsBetas       = zeros(this.nParams, this.NBETA);  
            this.paramsPopulations = zeros(this.nParams, this.NPOP);   
            this.paramsSigmas      = zeros(this.nParams, 1);       
            this.annealingAvpar    = zeros(this.nParams, 1);
            this.annealingSdpar    = zeros(this.nParams, 1);
            this.annealingInitz    = zeros(this.nParams, 1);
            this.bestFitParams     = zeros(this.nParams, 1);
            
            this.lpBetas       = zeros(this.NBETA, 1);
            this.lpPopulations = zeros(this.NPOP, 1);            
            
            
            this.paramsHist = zeros(this.nParams, this.NPOP*this.nProposalsQC); % for histogram of parameters
            this.logProbQC  = zeros(this.NPOP, this.nProposalsQC);              % qc on the lprob
            this.stdOfError = zeros(this.NPOP*this.nProposalsQC, 1);            % follow the standard dev of error            
            
            %% %%%%%%%%%%%%%%%%%%%
            %% initialize the MCMC
            %% %%%%%%%%%%%%%%%%%%%
            
            for m = 1:this.NPOP
                this.paramsSigmas = (this.paramsData.max - this.paramsData.min)/10.0;
                for k = 1:this.nParams
                    if (this.paramsData.fixed(k))
                        this.paramsPopulations(k,m) = this.paramsData.fixedValue(k);
                    else
                        this.paramsPopulations(k,m) = this.paramsData.mean(k) + this.paramsData.std(k)*randn(1,1);
                        while (this.paramsPopulations(k,m)<this.paramsData.min(k) || this.paramsPopulations(k,m)>this.paramsData.max(k))
                            this.paramsPopulations(k,m) = this.paramsData.mean(k) + this.paramsData.std(k)*randn(1,1);
                        end
                    end
                end
            end
        end  
        function [parmax,avpar,this] = runMcmc(this)
            %% MCMC (Markov Chain Monte-Carlo) has the machinery to do a simple Bayesian estimation
            %
            %  after J. S. Shimony, Mar, 2014

            %% %%%%%%%%%%%%%%%%%%%%%%%%
            %% annealing /burn in phase
            %% %%%%%%%%%%%%%%%%%%%%%%%%

            lp0 = nan;
            for b = 1:this.NBETA  
                
                beta_ = (1.0/(this.NBETA-1.0))*b; 
                this.printBeta(b, beta_, lp0);
                parn = this.annealingInitz;
                this.lpBetas(b) = 0.0;

                %% population loop
                for m = 1:this.NPOP
                    
                    ptmp = this.paramsPopulations(:,m);
                    lp0 = this.logProbability(ptmp, beta_, 1);
                    this.lpPopulations(m) = lp0;

                    for j = 1:this.NANNEAL
                        for k = 1:this.nParams
                            
                            if (this.paramsData.fixed(k)) 
                                continue; 
                            else
                                ptmp = this.paramsPopulations(:,m);
                                nprop = 0;
                                while (nprop < 50)
                                    nprop = nprop + 1;
                                    dpar = this.paramsSigmas(k)*randn(1,1);
                                    if (ptmp(k)+dpar>=this.paramsData.min(k) && ptmp(k)+dpar<=this.paramsData.max(k)) 
                                        break; 
                                    end
                                end
                                if (nprop < 50)
                                    ptmp(k) = ptmp(k) + dpar;
                                else
                                    fprintf('Warning: nprop too large, par %d val %f sigma %f\n',k,ptmp(k),this.paramsSigmas(k));
                                end
                                % new lp
                                lp1 = this.logProbability(ptmp, beta_, 1);

                                % metropolis acceptance
                                dele = lp1 - lp0;
                                if (dele > 0.0 || rand(1,1) < exp(dele))
                                    lp0 = lp1;
                                    this.lpPopulations(m) = lp1;
                                    parn(k) = parn(k) + 1;
                                    this.paramsPopulations(:,m) = ptmp;
                                end
                            end % end if (this.paramsData.fixed(k))
                        end % end for k = 1:this.nParams
                    end % end for j = 1:this.NANNEAL
                end % end for m = 1:this.NPOP

                this = this.gatherAnnealingStats(b);
                this = this.BretthorstAdjustments(parn);    
                this.printAnnealing(b, parn);
                this = this.replacePoorMembers;

            end % end for b = 1:this.NBETA 

            if (this.showPlots); this.plotAnnealing; end
            
            %% %%%%%%%%%%%%%%%%%%%%%%%
            %% proposal/sampling phase
            %% %%%%%%%%%%%%%%%%%%%%%%%

            for m = 1:this.NPOP
                
                ptmp = this.paramsPopulations(:,m);
                this.bestFitParams = this.paramsPopulations(:,m);
                lp0 = this.logProbability(ptmp, 1.0, 1);
                this.lpPopulations(m) = lp0;
                lpmax = lp0;
                parn = this.annealingInitz; 
                
                for j = 1:this.NPROPOSALS
                    for k = 1:this.nParams
                        
                        if (this.paramsData.fixed(k)) 
                            continue; 
                        else
                            ptmp = this.paramsPopulations(:,m);
                            nprop = 0;
                            while (nprop < 50)
                                nprop = nprop + 1;
                                dpar = this.paramsSigmas(k)*randn(1,1);
                                if (ptmp(k)+dpar>=this.paramsData.min(k) && ptmp(k)+dpar<=this.paramsData.max(k)) 
                                    break; 
                                end
                            end
                            if (nprop < 50)
                                ptmp(k) = ptmp(k) + dpar;
                            else
                                fprintf('Warning: nprop too large, par %d val %f sigma %f\n',k,ptmp(k),this.paramsSigmas(k));
                            end
                            
                            % new lp
                            lp1 = this.logProbability(ptmp, 1.0, 1);

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
                        end % end if (this.paramsData.fixed(k))
                    end % end for k = 1:this.nParams
                    
                    [lp1, this] = this.updateHistograms(j, m, lp0, lp1, ptmp);
                    
                end % end for j = 1:this.NPROPOSALS
            end % end for m = 1:this.NPOP
            
            %% %%%%%%%%%
            %% reporting
            %% %%%%%%%%%
            
            parmax = this.bestFitParams;
            avpar  = this.annealingAvpar;
            this.lpFinal = lp1;

            if (this.showBestFit);    this.printBestFit; end
            if (this.showFinalStats); this.printFinalStats; end
            if (this.showPlots);      this.histParametersDistributions; end
            if (this.showPlots);      this.histStdOfError; end
            if (this.showPlots);      this.plotLogProbabilityQC; end
        end    
        function lprob               = logProbability(this, paramsVec, beta_, lpFlag)
            %% LOGPROBABILITY
            %
            %  lpFlag: -1  return Q
            %           0  return log of normalized Q
            %           1  return with beta and priors
            %
            %  from J. S. Shimony, Mar, 2014
            
            lprob = this.bayesianProblem_.sumSquaredErrors(paramsVec);

            if (lpFlag == -1)
                return
            end

            % use for sigma fit
            if (paramsVec(1) < 0.0)
                lprob = lprob / (-2.0*paramsVec(1)^2);
                lprob = lprob + (-0.5*this.nSamples)*log(2.0*pi*paramsVec(1)^2);
            else  % no sigma, use t distribution Jeffrey Prior
                lprob = -0.5*this.nSamples*log(0.5*lprob);
            end

            % add in beta and pretest probabilities 
            if (lpFlag == 1)
                % beta only operates on likelihoods
                args  = this.PARPEN - (paramsVec - this.paramsData.mean).^2 ./ (2 * this.paramsData.std.^2);
                args  = double(~this.paramsData.fixed) .* args;
                lprob = beta_ * lprob + sum(args);
            end
        end

        function printBestFit(this)
            
            % print best fit member of population
            for k = 1:this.nParams
                fprintf('BEST-FIT    param %s value   %f\n', this.paramIndexToLabel(k), this.bestFitParams(k));
            end
        end
        function printFinalStats(this)
            
            avpar = this.annealingInitz;
            sdpar = this.annealingInitz;
            for m = 1:this.NPOP
                ptmp = this.paramsPopulations(:,m);
                for k = 1:this.nParams
                    avpar(k) = avpar(k) + ptmp(k);
                    sdpar(k) = sdpar(k) + ptmp(k)^2;
                end
            end
            for k = 1:this.nParams
                avpar(k) = avpar(k)/this.NPOP;
                sdpar(k) = (sdpar(k) - this.NPOP*avpar(k)^2)/(this.NPOP - 1);
                if (sdpar(k) > 0.0)
                    sdpar(k) = sqrt(sdpar(k));
                else
                    sdpar(k) = 0.0;
                end
            end
            
            avpar = mean(this.paramsPopulations'); %#ok<UDIM>
            sdpar = std(this.paramsPopulations'); %#ok<UDIM>
            for k = 1:this.nParams
                fprintf('FINAL STATS param %s average %f std %f\n', this.paramIndexToLabel(k), avpar(k), sdpar(k));
            end
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
        function histStdOfError(this)
            figure;
            hist(this.stdOfError, this.NBINS);
            xlabel('std. dev. of error');
        end
        function plotAnnealing(this)
            
            % plotting on annealing phase
            % for k = 1:NPAR
            %     figure;
            %     plot(1:NBETA, this.paramsBetas(k,:));
            %     title(['Parameter ',num2str(k),' vs temperature']);
            % end
            
            figure;
            plot(1:this.NBETA, this.lpBetas);
            xlabel('beta (1/temp)');
            ylabel('log(prob)');
        end
        function plotLogProbabilityQC(this)
            figure;
            hold on;
            for k = 1:this.NPOP
                plot(1:this.nProposalsQC, this.logProbQC(k,:), 'Color', this.colorVariation(k, this.NPOP));
            end
            hold off;
            title('population index:  red->blue');
            xlabel('proposal#');
            ylabel('log(prob)')
        end
    end
    
    %% PRIVATE
    
    properties (Access = 'private')
        bayesianProblem_
    end
    
    methods (Access = 'private')
        function              printBeta(this, b, beta, lp0)
            if (this.verbosity < eps)
                return; end
            fprintf('annealing step %d beta (1/temp) %d logProb0 %d\n', b, beta, lp0); 
        end        
        function              printAnnealing(this, b, parn)
            if (this.verbosity < 0.2)
                return; end
            if (1                        == b || ...
                this.FRACPEEK*this.NBETA == b || ...
                this.NBETA               == b)
                for k = 1:this.nParams
                    fprintf('annealing step %d par %s avg %f std %f sigma %f acc %f\n',...
                        b,this.paramIndexToLabel(k),this.annealingAvpar(k),this.annealingSdpar(k),this.paramsSigmas(k),parn(k)/(this.NANNEAL*this.NPOP));
                end
            end
        end
        
        function [this,lp1] = gatherAnnealingStats(this, b)
            this.annealingAvpar = this.annealingInitz;
            this.annealingSdpar = this.annealingInitz;
            for m = 1:this.NPOP
                ptmp = this.paramsPopulations(:,m);
                for k = 1:this.nParams
                    this.annealingAvpar(k) = this.annealingAvpar(k) + ptmp(k);
                    this.annealingSdpar(k) = this.annealingSdpar(k) + ptmp(k)^2;
                end
                beta_ = (1.0/(this.NBETA-1.0))*b; 
                lp1 = this.logProbability(ptmp, beta_, 0);
                this.lpBetas(b) = this.lpBetas(b) + lp1/this.NPOP;
            end
            for k = 1:this.nParams
                this.annealingAvpar(k) = this.annealingAvpar(k)/this.NPOP;
                this.annealingSdpar(k) = (this.annealingSdpar(k) - this.NPOP*this.annealingAvpar(k)^2)/(this.NPOP - 1);
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
                if (parn(k) < 0.1*this.NANNEAL*this.NPOP)
                        this.paramsSigmas(k) = 0.1*this.paramsSigmas(k);
                elseif (parn(k) < 0.2*this.NANNEAL*this.NPOP)
                    if (this.paramsSigmas(k) > 0.1*this.annealingSdpar(k))
                        this.paramsSigmas(k) = 0.8*this.paramsSigmas(k);
                    end
                elseif (parn(k) > 0.8*this.NANNEAL*this.NPOP)
                    if (this.paramsSigmas(k) < 5.0*this.annealingSdpar(k))
                        this.paramsSigmas(k) = 10.0*this.paramsSigmas(k);
                    end
                elseif (parn(k) > 0.55*this.NANNEAL*this.NPOP)
                    if (this.paramsSigmas(k) < 5.0*this.annealingSdpar(k))
                        this.paramsSigmas(k) = 2.0*this.paramsSigmas(k);
                    end
                elseif (parn(k) > 0.3*this.NANNEAL*this.NPOP)
                    if (this.paramsSigmas(k) < 5.0*this.annealingSdpar(k))
                        this.paramsSigmas(k) = 1.1*this.paramsSigmas(k);
                    end
                end
                if (this.paramsSigmas(k) > (this.paramsData.max(k) - this.paramsData.min(k)))
                    this.paramsSigmas(k) =  this.paramsData.max(k) - this.paramsData.min(k);
                end           
            end  
        end
        function  this      = replacePoorMembers(this)
            for k = 1:this.NPOPREP
                lpmax = -this.LRG;
                lpmin = this.LRG;
                ilpmax = 0;
                ilpmin = 0;
                for m = 1:this.NPOP
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
            map  = this.paramsData.paramsMap;
            keys = map.keys;
            lbl  = keys{idx};
            if (isnumeric(lbl))
                lbl = num2str(lbl);
            end
        end
        function [lp1,this] = updateHistograms(this, j, m, lp0, lp1, ptmp)
            if (mod(j,1/this.FRACPEEK) == 0)
                this.logProbQC(m, j*this.FRACPEEK) = lp0;
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

