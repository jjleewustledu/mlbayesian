classdef BlomqvistMcmcProblem < mlbayesian.AbstractMcmcProblem
	%% BLOMQVISTMCMCPROBLEM is for testing AbstractMcmcProblem implementations

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties
        showPlots = true
        baseTitle = 'BlomqvistMcmcProblem'
        xLabel    = 'times'
        yLabel    = 'wellcounts'
        
        Ca
        VB = 0.0369
        K04 = 0.26057
        dt = 1
    end
    
	methods 
        function this = BlomqvistMcmcProblem(t, y, ca)
            this = this@mlbayesian.AbstractMcmcProblem(t, y);
            this.Ca = ca;
        end
        function this = estimateParameters(this)
            %% ESTIMATEPARAMETERS
            %  expected_parameters = [0.2606 4.7967e-04 0.0348 0.0025 4.2528e-04 0]; % k04 k12 k21 k32 k43 t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('k04') = struct('fixed', 1, 'min', 1e-3, 'mean', this.K04,   'max', 10); 
            map('k12') = struct('fixed', 1, 'min', 1e-3, 'mean', 0.00047967, 'max', 0.05);
            map('k21') = struct('fixed', 0, 'min', 1e-2, 'mean', 0.034772,   'max', 0.1);
            map('k32') = struct('fixed', 1, 'min', 1e-3, 'mean', 0.0025173,  'max', 10);
            map('k43') = struct('fixed', 1, 'min', 1e-4, 'mean', 0.00042528, 'max', 0.005);  
            map('t0' ) = struct('fixed', 1, 'min',    0, 'mean', 0,          'max', 10);  

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('k04'), this.finalParams('k12'), this.finalParams('k21'), ...
                this.finalParams('k32'), this.finalParams('k43'), this.finalParams('t0'));
        end
        function ed   = estimateDataFast(this, k04, k12, k21, k32, k43, t0)
        end   
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

