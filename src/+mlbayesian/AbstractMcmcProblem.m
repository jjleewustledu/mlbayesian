classdef AbstractMcmcProblem < mlbayesian.AbstractBayesianProblem 
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
 	 

    properties (Abstract)
        baseTitle
        xLabel
        yLabel
    end
    
    properties (Dependent)
        length % of dependent_data = f(time_interpolants), which must have the same array sizes
        timeInterpolants
    end
    
    methods %% GET
        function le = get.length(this)
            le = length(this.independentData);
        end
        function ti = get.timeInterpolants(this)
            ti = this.independentData;
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
        function        plotOri(this)
            figure
            plot(this.timeInterpolants, this.dependentData, 'k', 'LineWidth', 2)
            title(this.baseTitle)
            xlabel(this.xLabel)
            ylabel(this.yLabel)
        end
        function        plotEstimate(this)
            figure
            plot(this.timeInterpolants, this.dependentData, 'k', ...
                 this.timeInterpolants, this.estimateData,  'k:', 'LineWidth', 2);
            title(sprintf('%s and Bayesian estimate', this.baseTitle));
            xlabel(this.xLabel)
            ylabel(this.yLabel)
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

