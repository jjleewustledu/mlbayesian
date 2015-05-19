classdef (Abstract) IPerfusionProblem  
	%% IPERFUSIONPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.5.0.197613 (R2015a) 
 	%  $Id$  	 

	properties (Abstract)
        baseTitle
        detailedTitle
        concentration_a
        concentration_obs
        mtt_a
        mtt_obs
        map 
        xLabel
        yLabel
    end

	methods (Static, Abstract)
        this  = load
        ci    = concentration_i
        args  = interpolateData
        [t,y] = shiftData
        [t,y] = shiftDataLeft
        [t,y] = shiftDataRight
        this  = simulateMcmc
    end 
    
    methods (Abstract)
        this = estimateAll(this)
        this = estimateParameters(this)
        ed   = estimateData(this)
        ed   = estimateDataFast(this)
        ci   = itsConcentration_i(this)
        this = save(this)
        this = saveas(this, fn)
        this = simulateItsMcmc(this)
        
               plotInitialData(this)
               plotProduct(this)
               plotParVars(this)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

