classdef KineticsMcmcProblem < mlbayesian.AbstractMcmcProblem
	%% KINETICSMCMCPROBLEM is for testing AbstractMcmcProblem implementations

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties
        showPlots = true
        baseTitle = 'KineticsMcmcProblem'
        xLabel    = 'times'
        yLabel    = 'wellcounts'
        
        Ca
        VB = 0.0369
        K04 = 0.26057
        dt = 1
    end
    
	methods 
        function this = KineticsMcmcProblem(t, y, ca)
            this = this@mlbayesian.AbstractMcmcProblem(t, y);
            this.Ca = ca;
        end
        function this = estimateParameters(this)
            %% ESTIMATEPARAMETERS
            %  expected_parameters = [0.2606 4.7967e-04 0.0348 0.0025 4.2528e-04 0]; % k04 k12 k21 k32 k43 t0
            
            import mlbayesian.*;
            map = containers.Map;
            map('k04') = struct('fixed', 1, 'min', 1e-3, 'mean', this.K04,   'max', 10); 
            map('k12') = struct('fixed', 0, 'min', 1e-4, 'mean', 0.00047967, 'max', 0.01);
            map('k21') = struct('fixed', 0, 'min', 1e-2, 'mean', 0.034772,     'max', 0.1);
            map('k32') = struct('fixed', 0, 'min', 1e-4, 'mean', 0.0025173,     'max', 10);
            map('k43') = struct('fixed', 0, 'min', 1e-4, 'mean', 0.00042528, 'max', 0.001);  
            map('t0' ) = struct('fixed', 1, 'min',    0, 'mean', 0,          'max', 10);  

            this.paramsManager = BayesianParameters(map);
%             this.paramsManager.NPROPOSALS = 200;
%             this.paramsManager.NPOP       = 100;
%             this.paramsManager.NPOPREP    =  10;
%             this.paramsManager.NBETA      = 100;
%             this.paramsManager.NANNEAL    =  20;
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;
        end
        function ed   = estimateData(this)
            ed = this.estimateDataFast( ...
                this.finalParams('k04'), this.finalParams('k12'), this.finalParams('k21'), ...
                this.finalParams('k32'), this.finalParams('k43'), this.finalParams('t0'));
        end
        function ed   = estimateDataFast(this, k04, k12, k21, k32, k43, t0)            
            k22 = k12 + k32;
            t = this.timeInterpolants;
            
            q2_ = this.VB * k21 * exp(-k22*t);
            q3_ = this.VB * k21 * k32 * (k22 - k43)^-1 * (exp(-k43*t) - exp(-k22*t));
            q4_ = this.VB * k21 * k32 * k43 * ( ...
                     exp(-k22*t)/((k04 - k22)*(k43 - k22)) + ...
                     exp(-k43*t)/((k22 - k43)*(k04 - k43)) + ...
                     exp(-k04*t)/((k22 - k04)*(k43 - k04)));
            q234 = conv(q2_ + q3_ + q4_, this.Ca);
            q234 = q234(1:length(t));
            Q0   = this.VB * this.Ca + q234;            
            
            t0_idx         = floor(t0/this.dt) + 1;
            ed             = zeros(size(Q0));
            ed(t0_idx:end) = Q0(1:end-t0_idx+1);
        end   
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

