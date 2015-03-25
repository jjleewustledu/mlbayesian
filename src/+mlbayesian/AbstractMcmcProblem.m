classdef AbstractMcmcProblem < mlbayesian.AbstractBayesianProblem & mlbayesian.IMcmcProblem
	%% ABSTRACTMCMCPROBLEM   
    %  Yet abstract:
    %      properties baseTitle, xLabel, yLabel, showPlots
    %      methods estimateParameters, estimateData, estimateDataFast

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 
    properties (Dependent)
        length % of dependent_data = f(time_interpolants), which must have the same array sizes
        timeInterpolants
        timeFinal
        
        NPROPOSALS % number of loops in parameter prob phase
        NPOP       % number of population
        NPOPREP    % number of population to replace
        NBETA      % number of temperature steps
        NANNEAL    % number of loops per annealing temp        
        
        annealingAvpar
        annealingSdpar
        annealingInitz
    end
    
    methods %% GET
        function le = get.length(this)
            le = length(this.independentData);
        end
        function ti = get.timeInterpolants(this)
            ti = this.independentData;
        end
        function tf = get.timeFinal(this)
            tf = this.independentData(end);
        end
        function n  = get.NPROPOSALS(this)
            n = this.mcmc.NPROPOSALS;
        end
        function n  = get.NPOP(this)
            n = this.mcmc.NPOP;
        end
        function n  = get.NPOPREP(this)
            n = this.mcmc.NPOPREP;
        end
        function n  = get.NBETA(this)
            n = this.mcmc.NBETA;
        end
        function n  = get.NANNEAL(this)
            n = this.mcmc.NANNEAL;
        end
        function a  = get.annealingAvpar(this)
            a = this.mcmc.annealingAvpar;
        end
        function a  = get.annealingSdpar(this)
            a = this.mcmc.annealingSdpar;
        end
        function a  = get.annealingInitz(this)
            a = this.mcmc.annealingInitz;
        end
    end

	methods 		  
 		function this = AbstractMcmcProblem(varargin) 
 			%% ABSTRACTMCMCPROBLEM 
 			%  Usage:  this = AbstractMcmcProblem() 

 			this = this@mlbayesian.AbstractBayesianProblem(varargin{:}); 
 		end 
        function this = runMcmc(this, paramsMap)
            %% RUNMCMC should be run from within method estimateParameters, implemented as below
            %  Usage:  map = containers.Map
            %          map('some_param') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = this.runMcmc(map)
            
            import mlbayesian.*;
            this.paramsManager = BayesianParameters(paramsMap);            
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;            
            
            if (this.showPlots)
                this.plotOri;
                this.plotEstimate;
            end
        end   
        function lp   = logProbability(this, paramsVec, beta, logProbabilityFlag)
            lp = this.mcmc.logProbability(paramsVec, beta, logProbabilityFlag);
        end
        function        printBestFit(this)
            this.mcmc.printBestFit;
        end
        function        printFinalStats(this)
            this.mcmc.printFinalStats;
        end
        function        histParametersDistributions(this) 
            this.mcmc.histParametersDistributions;
        end
        function        histStdOfError(this) 
            this.mcmc.histStdOfError;
        end
        function        plotAnnealing(this) 
            this.mcmc.plotAnnealing;
        end
        function        plotLogProbabilityQC(this) 
            this.mcmc.plotLogProbabilityQC;
        end    
        function        plotOri(this)
            figure
            plot(this.timeInterpolants, this.dependentData, 'o')
            title(this.baseTitle)
            xlabel(this.xLabel)
            ylabel(this.yLabel)
        end
        function        plotEstimate(this)
            figure
            plot(this.timeInterpolants, this.dependentData, 'o', ...
                 this.timeInterpolants, this.estimateData,  '-');
            title(sprintf('%s and Bayesian estimate', this.baseTitle));
            xlabel(this.xLabel)
            ylabel(this.yLabel)
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

