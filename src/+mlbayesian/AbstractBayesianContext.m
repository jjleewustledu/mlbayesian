classdef (Abstract) AbstractBayesianContext
	%% ABSTRACTBAYESIANCONTEXT is the context for AbstractBayesianStrategy 

	%  $Revision$
 	%  was created 23-Nov-2015 17:06:10
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64. 	
    
    properties (Dependent)
        theStrategy
    end
    
    methods %% GET
        function ts = get.theStrategy(this)
            ts = this.theStrategy_;
        end
    end
    
    %% PRIVATE
    
    properties (Access = 'private')
        theStrategy_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

