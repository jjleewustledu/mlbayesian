classdef McmcKernel 
	%% MCMCKERNEL is designed for recoding for compiled objects.  It traces its legacy to Joshua S. Shimony, 
    %  G. Larry Bretthorst and Edwin Jaynes.

	%  $Revision$
 	%  was created 14-Dec-2017 18:02:49 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John J. Lee and Joshua S. Shimony.
 	
    properties (Constant)
        FRAC_POPREP =  0.1
        FRACPEEK    =  0.2
        MAX_PROP    = 50
        N_SIGMAS    = 10   % of normal distribution to use for range [min, max] of parameters
    end
    
	properties 		
        annealingAvpar
        annealingSdpar
        annealingInitz
        bestFitParams
        isfinished
        meanParams
        paramsBetas
        paramsHist
        paramsPopulations
        priorSigmas
        showAnnealing
        showBeta
        stdOfError
        stdParams
        verbosity
        
        lpQC        
        lpBetas
        lpPopulations
        lpFinal
    end
    
    methods (Static)
        function this = main(varargin)
            this = mlbayesian.McmcKernel(varargin{:});
            this = this.runMcmc;
        end
    end

	methods
 		function this = McmcKernel(mcmcs)
 			%% MCMCKERNEL
            %  @param mcmcStruct is a struct
            
            assert(isstruct(mcmcs));
            this.mcmcStruct_    = mcmcs;   
            this.paramsPenalty_ = this.mcmcStruct_.paramsPenalty;         
            this.paramsStruct_  = this.mcmcStruct_.solverParameters; % another struct, cached for speed
            
            % problems here affect replacePoorMembers, logProbability
            assert(all(size(this.paramsStruct_.mean_)      == [this.nParams 1])); 
            assert(all(size(this.paramsStruct_.std_)       == [this.nParams 1]));
            assert(all(size(this.paramsStruct_.fixed)      == [this.nParams 1]));
            assert(all(size(this.paramsStruct_.fixedValue) == [this.nParams 1]));
            
            this.annealingAvpar    = zeros(this.nParams, 1);
            this.annealingSdpar    = zeros(this.nParams, 1);
            this.annealingInitz    = zeros(this.nParams, 1);
            this.bestFitParams     = zeros(this.nParams, 1);  
            this.meanParams        = zeros(this.nParams, 1);
            this.paramsBetas       = zeros(this.nParams, this.nBeta);  
            this.paramsHist        = zeros(this.nParams, this.nPop*this.nProposalsQC); % for histogram of parameters  
            this.paramsPopulations = zeros(this.nParams, this.nPop);   
            this.priorSigmas       = zeros(this.nParams, 1);            
            this.showAnnealing     = this.mcmcStruct_.showAnnealing;
            this.showBeta          = this.mcmcStruct_.showBeta;
            this.stdOfError        = zeros(this.nPop*this.nProposalsQC, 1); % follow the standard deviation of error 
            this.stdParams         = zeros(this.nParams, 1);
            this.verbosity         = this.mcmcStruct_.verbosity;    
            
            this.lpQC          = zeros(this.nPop, this.nProposalsQC); % qc on the lprob   
            this.lpBetas       = zeros(this.nBeta, 1);
            this.lpPopulations = zeros(this.nPop, 1);            
            
            
            %% initialize the Mcmc
            if (this.verbosity > eps)
                fprintf('mlbayesian.McmcKernel.ctor:  initializing McmcKernel'); end
            this.priorSigmas = (this.paramsStruct_.max_ - this.paramsStruct_.min_)/this.N_SIGMAS;
            for m = 1:this.nPop
                for k = 1:this.nParams
                    if (this.paramsStruct_.fixed(k))
                        this.paramsPopulations(k,m) = this.paramsStruct_.fixedValue(k);
                    else
                        this.paramsPopulations(k,m) = this.paramsStruct_.mean_(k) + this.paramsStruct_.std_(k)*randn(1,1);
                        while (this.paramsPopulations(k,m)<this.paramsStruct_.min_(k) || ...
                               this.paramsPopulations(k,m)>this.paramsStruct_.max_(k))
                            this.paramsPopulations(k,m) = this.paramsStruct_.mean_(k) + this.paramsStruct_.std_(k)*randn(1,1);
                        end
                    end
                end
            end 
            this.isfinished = false;
        end
        
        function n = nAnneal(this)
            n = this.paramsStruct_.nAnneal;
        end
        function n = nBeta(this)
            n = this.paramsStruct_.nBeta;
        end
        function n = nParams(this)
            n = length(this.paramsStruct_.mean_);
        end
        function n = nPop(this)
            n = this.paramsStruct_.nPop;
        end
        function n = nPopRep(this)
            n = this.FRAC_POPREP * this.nPop;
        end
        function n = nProposals(this)
            n = this.paramsStruct_.nProposals;
        end
        function n = nProposalsQC(this)
            n = ceil(this.FRACPEEK * this.nProposals);
        end  
        function n = nSamples(this)
            n = this.paramsStruct_.nSamples;
        end        
        function this = runMcmc(this)
            %% RUNMCMC (Markov Chain Monte-Carlo) has the machinery to do a simple Bayesian estimation
            %
            %  after J. S. Shimony, Mar, 2014

            %% %%%%%%%%%%%%%%%%%%%%%%%
            %% annealing/burn-in phase
            %% %%%%%%%%%%%%%%%%%%%%%%%
            
            if (this.verbosity > eps)
                fprintf('mlbayesian.McmcKernel.runMcmc:  annealing/burn-in'); 
            end
            dpar = nan;
            lp0  = nan;
            lp1  = nan;
            for b = 1:this.nBeta  
                
                beta_ = (1/(this.nBeta - 1))*b; 
                this.printBeta(b, beta_, lp0);
                parn = this.annealingInitz;
                this.lpBetas(b) = 0;

                %% population loop
                
                for m = 1:this.nPop
                    
                    ptmp = this.paramsPopulations(:,m);
                    lp0  = this.logProbability(ptmp, beta_, 1);
                    this.lpPopulations(m) = lp0;

                    for j = 1:this.nAnneal
                        for k = 1:this.nParams
                            
                            if (this.paramsStruct_.fixed(k)) 
                                continue; 
                            else
                                ptmp = this.paramsPopulations(:,m);
                                nprop = 0;
                                while (nprop < this.MAX_PROP)
                                    nprop = nprop + 1;
                                    dpar = this.priorSigmas(k)*randn(1,1);
                                    if (ptmp(k)+dpar >= this.paramsStruct_.min_(k) && ptmp(k)+dpar <= this.paramsStruct_.max_(k)) 
                                        break; 
                                    end
                                end
                                if (nprop < this.MAX_PROP)
                                    ptmp(k) = ptmp(k) + dpar;
                                else
                                    if (1 == this.verbosity)
                                        fprintf('Warning: nprop %i, par %g, val %g, sigma %g\n', ...
                                            int64(nprop), k, ptmp(k), this.priorSigmas(k)); 
                                    end
                                end
                                % new lp
                                [lp1,ptmp] = this.logProbability(ptmp, beta_, 1);

                                % metropolis acceptance
                                dele = lp1 - lp0;
                                if (dele > 0 || rand(1,1) < exp(dele))
                                    lp0 = lp1;
                                    this.lpPopulations(m) = lp1;
                                    parn(k) = parn(k) + 1;
                                    this.paramsPopulations(:,m) = ptmp;
                                end
                            end % end if (this.paramsStruct_.fixed(k))
                        end % end for k = 1:this.nParams
                    end % end for j = 1:this.nAnneal
                end % end for m = 1:this.nPop
                
                %%

                this = this.gatherAnnealingStats(b);
                this = this.BretthorstAdjustments(parn);    
                this.printAnnealing(b, parn);
                this = this.replacePoorMembers;

            end % end for b = 1:this.nBeta 
            
            %% %%%%%%%%%%%%%%%%%%%%%%%
            %% proposal/sampling phase
            %% %%%%%%%%%%%%%%%%%%%%%%%

            if (this.verbosity > eps)
                fprintf('mlbayesian.McmcKernel.runMcmc:  proposal/sampling'); 
            end
            for m = 1:this.nPop
                
                ptmp = this.paramsPopulations(:,m);
                this.bestFitParams = this.paramsPopulations(:,m);
                [lp0,ptmp] = this.logProbability(ptmp, 1, 1);
                this.lpPopulations(m) = lp0;
                lpmax = lp0;
                parn = this.annealingInitz; 
                
                %% proposal (importance sampling) loop
                
                for j = 1:this.nProposals
                    for k = 1:this.nParams
                        
                        if (this.paramsStruct_.fixed(k)) 
                            continue; 
                        else
                            ptmp = this.paramsPopulations(:,m);
                            nprop = 0;
                            while (nprop < this.MAX_PROP)
                                nprop = nprop + 1;
                                dpar = this.priorSigmas(k)*randn(1,1);
                                if (ptmp(k)+dpar>=this.paramsStruct_.min_(k) && ptmp(k)+dpar<=this.paramsStruct_.max_(k)) 
                                    break; 
                                end
                            end
                            if (nprop < this.MAX_PROP)
                                ptmp(k) = ptmp(k) + dpar;
                            else
                                if (1 == this.verbosity)
                                    fprintf('Warning: nprop %i, par %g, val %g, sigma %g\n', ...
                                        int64(nprop),k,ptmp(k),this.priorSigmas(k)); 
                                end
                            end
                            
                            % new lp
                            [lp1,ptmp] = this.logProbability(ptmp, 1, 1);

                            % store best for future use
                            if (lp1 > lpmax)
                                lpmax = lp1;
                                this.bestFitParams = ptmp;
                            end

                            % metropolis acceptance
                            dele = lp1 - lp0;
                            if (dele > 0 || rand(1,1) < exp(dele))
                                lp0 = lp1;
                                this.lpPopulations(m) = lp1;
                                parn(k) = parn(k) + 1;
                                this.paramsPopulations(:,m) = ptmp;
                            end
                        end % end if (this.paramsStruct_.fixed(k))
                    end % end for k = 1:this.nParams
                    
                    [lp1, this] = this.updateHistograms(j, m, lp0, lp1, ptmp);
                    
                end % end for j = 1:this.nProposals
                
                %%
                
            end % end for m = 1:this.nPop
            
            %%
            
            this.lpFinal = lp1;
            this = this.finalize;
        end
        function this = finalize(this)
            this.isfinished = true;
            
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
            for k = 1:this.nParams
                this.meanParams(k) = avpar(k);
                this.stdParams(k)  = sdpar(k);
            end
        end
        function [lprob,paramsVec] = logProbability(this, paramsVec, beta_, lpFlag)
            %% LOGPROBABILITY
            %
            %  lpFlag: -1  return Q
            %           0  return log of normalized Q
            %           1  return with beta and priors
            %
            %  from J. S. Shimony, Mar, 2014
            
            paramsVec = this.mcmcStruct_.adjustParams(paramsVec);
            lprob     = this.mcmcStruct_.objectiveFunc(paramsVec);
            
            if (lpFlag == -1)
                return
            end

            % use log t-distribution for parameter degrees of freedom; cf. Lee MRM 2010 eqn 10
            lprob = -0.5*this.nSamples*log(0.5*lprob);

            % add in beta and pretest probabilities
            if (lpFlag == 1)
                % beta only operates on likelihoods
%                 m_    = this.paramsStruct_.mean_;
%                 s_    = this.paramsStruct_.std_;
%                 pv_   = paramsVec;
%                 args  = this.paramsPenalty_ - (pv_ - m_).^2 ./ (2 * s_.^2);
%                 args  = args(~this.paramsStruct_.fixed);
%                 lprob = beta_*lprob + sum(args);
                
                len   = length(paramsVec);
                lprob = beta_ * lprob;
                for f = 1:len
                   if (~this.paramsStruct_.fixed(f))
                       lprob = lprob + this.paramsPenalty_ - (paramsVec(f) - this.paramsStruct_.mean_(f))^2/(2.0*this.paramsStruct_.std_(f)^2);
                   end
                end
            end
        end
        function fprintf(this)
            fprintf(['mlbayesian.McmcKernel:\n' this.sprintf]);
        end
        function s = sprintf(this)
            s = '';
            for k = 1:this.nParams
                s = [s sprintf('FINAL STATS param %s \tmode %g \tmean %g \tstd %g\n', ...
                     this.mcmcStruct_.parameterIndexToName(k), ...
                     this.bestFitParams(k), this.meanParams(k), this.stdParams(k))]; %#ok<AGROW>
            end
        end
    end 
    
    %% PRIVATE
    
    properties (Access = 'private')
        mcmcStruct_
        paramsStruct_
        paramsPenalty_
    end
    
    methods (Access = 'private')
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
                if (this.paramsStruct_.fixed(k)) 
                    continue;
                end
                
                if (parn(k) < 0.1*this.nAnneal*this.nPop)
                        this.priorSigmas(k) = 0.2 *this.priorSigmas(k);
                elseif (parn(k) < 0.2*this.nAnneal*this.nPop)
                    if (this.priorSigmas(k) > 0.05*this.annealingSdpar(k))
                        this.priorSigmas(k) =      this.priorSigmas(k)/1.2;
                    end
                elseif (parn(k) > 0.8*this.nAnneal*this.nPop)
                    if (this.priorSigmas(k) <     this.annealingSdpar(k))
                        this.priorSigmas(k) = 5.0 *this.priorSigmas(k);
                    end
                elseif (parn(k) > 0.6*this.nAnneal*this.nPop)
                    if (this.priorSigmas(k) <     this.annealingSdpar(k))
                        this.priorSigmas(k) = 2.0*this.priorSigmas(k);
                    end
                elseif (parn(k) > 0.3*this.nAnneal*this.nPop)
                    if (this.priorSigmas(k) <     this.annealingSdpar(k))
                        this.priorSigmas(k) = 1.1*this.priorSigmas(k);
                    end
                end
                if (this.priorSigmas(k) > (this.paramsStruct_.max_(k) - this.paramsStruct_.min_(k)))
                    %warning('mlbayesian:parameterOutOfBounds', ...
                    %    'McmcKernel.BretthorstAdjustments.priorSigmas(%i)->%g is too large\n', int64(k), this.priorSigmas(k));
                    this.priorSigmas(k) =  this.paramsStruct_.max_(k) - this.paramsStruct_.min_(k);
                end           
            end  
        end
        function  this      = replacePoorMembers(this)
            for k = 1:this.nPopRep
                lpmax = -Inf;
                lpmin =  Inf;
                ilpmax = 0;
                ilpmin = 0;
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
                
                %% assure that these are not repeated;
                %try
                    this.lpPopulations(ilpmin) = 0.0;
                    this.lpPopulations(ilpmax) = 0.0;
                    this.paramsPopulations(:,ilpmin) = this.paramsPopulations(:,ilpmax);
                %catch ME
                %    fprintf('\n');
                %    fprintf('    ---------\n');
                %    fprintf('    Error ID:\n');
                %    fprintf('    ---------\n');
                %    fprintf('    ''mlbayesian:possibleMissingData''\n\n');
                %    fprintf('    --------------\n');
                %    fprintf('    Error Details:\n');
                %    fprintf('    --------------\n');
                %    fprintf('    McmcKernel.replacePoorMembers.paramsStruc_.fixedValue -> %s\n', ...
                %        mat2str(this.paramsStruct_.fixedValue));
                %    handexcept(ME);
                %end
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
        function              printBeta(this, b, beta, lp0)
            if (~this.showBeta); return; end
            fprintf('annealing step %g beta (1/temp) %g logProb0 %g\n', b, beta, lp0); 
        end        
        function              printAnnealing(this, b, parn)
            if (~this.showAnnealing); return; end
            if ( ceil(this.FRACPEEK*this.nBeta) == b)
                fprintf('\n');
                for k = 1:this.nParams
                    fprintf('annealing step %g param %3s mean %g\t std %g\t sigma_{prior} %g\t acc %g\n',...
                        b, ...
                        this.mcmcStruct_.parameterIndexToName(k), ...
                        this.annealingAvpar(k), ...
                        this.annealingSdpar(k), ...
                        this.priorSigmas(k), ...
                        parn(k)/(this.nAnneal*this.nPop));
                end
                fprintf('\n');
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

