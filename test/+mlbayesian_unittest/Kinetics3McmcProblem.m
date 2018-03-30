classdef Kinetics3McmcProblem < mlbayesian.AbstractMcmcProblem
	%% KINETICS3MCMCPROBLEM is for testing AbstractMcmcProblem implementations

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties
        baseTitle = 'Kinetics3McmcProblem'
        xLabel    = 'times'
        yLabel    = 'wellcounts'
        
        Ca
        VB = 0.0369
        K04 = 0.26057
        dt = 1
    end
    
	methods 
        function this = Kinetics3McmcProblem(t, y, ca)
            this = this@mlbayesian.AbstractMcmcProblem(t, y);
            this.Ca = ca;
        end
        function this = estimateParameters(this)
            %% ESTIMATEPARAMETERS
            %  expected_parameters = [0.2606 4.7967e-04 0.0348 0.00042528 4.3924]; % k04 k12 k21 k32 t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('k04') = struct('fixed', 1, 'min', 0.01, 'mean', this.K04,    'max', 1); 
            map('k12') = struct('fixed', 0, 'min', 2e-3, 'mean', 0.0010001,   'max', 0.02);
            map('k21') = struct('fixed', 0, 'min', 0.04, 'mean', 0.0650655,   'max', 0.09);
            map('k32') = struct('fixed', 0, 'min', 1e-3, 'mean', 1.12944e-06, 'max', 8);
            map('t0' ) = struct('fixed', 0, 'min',    0, 'mean', 30,          'max', 30);  

            this.paramsManager = BayesianParameters(map);
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('k04'), this.finalParams('k12'), this.finalParams('k21'), ...
                this.finalParams('k32'),                          this.finalParams('t0'));
        end
        function Q    = estimateDataFast(this, k04, k12, k21, k32, t0)            
            k22 = k12 + k32;
            t = this.timeInterpolants;
            
            t0_idx               = floor(t0/this.dt) + 1;
            cart                 = this.Ca(end) * ones(1, this.length);
            cart(1:end-t0_idx+1) = this.Ca(t0_idx:end); 
            
            q2_ = this.VB * k21 * exp(-k22*t);
            q3_ = this.VB * k21 * k32 * (k22 - k04)^-1 * (exp(-k04*t) - exp(-k22*t));
            q23 = conv(q2_ + q3_, cart);
            q23 = q23(1:length(t));
            Q   = this.VB * cart + q23;
        end   
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

