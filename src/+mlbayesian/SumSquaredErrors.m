classdef (Abstract) SumSquaredErrors < handle
	%% SUMSQUAREDERRORS  

	%  $Revision$
 	%  was created 28-Aug-2018 17:07:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.

    methods (Static)
        function this = CreateSumSquaredErrors(client, QType)
            assert(ischar(QType));
            import mlbayesian.*;
            switch (QType)
                case 'SumSquaredResiduals'
                    this = SumSquaredResiduals(client);
                case 'SumSquaredWeightedResiduals'
                    this = SumSquaredWeightedResiduals(client);
                case 'SumSquaredJeffreysResiduals'
                    this = SumSquaredJeffreysResiduals(client);
            end
        end
    end
    
    %% PROTECTED
    
	methods (Access = protected)
 		function this = SumSquaredErrors(varargin)
 			%% SUMSQUAREDERRORS
 			%  @param .
 	
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'client', @(x) isa(x, 'mlbayesian.IMcmcStrategy'));
            parse(ip, varargin{:});
            this.client_ = ip.Results.client;
 		end
    end 
    
    properties (Access = protected)
        client_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

