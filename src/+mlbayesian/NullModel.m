classdef NullModel < mlbayesian.IBretthorstModel & mlio.AbstractIO
	%% NULLMODEL  

	%  $Revision$
 	%  was created 26-Dec-2017 23:28:05 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee. 	
    
    properties (Dependent)
        kernel
 		independentData
        dependentData
    end
    
	properties
        classOfSolver
        useSynthetic          
        %M        = mlbayesian.McmcKernel.N_SIGMAS/2; % match the prior width to the annealing width
        M          = 50; % noninformative prior width
        nAnneal    = 20
        nBeta      = 50
        nProposals = 100
    end

    methods (Static)
        function name = parameterIndexToName(varargin)
            name = '';
        end
    end
    
	methods
        
        %% GET, SET
        
        function g    = get.kernel(this)
            g = this.kernel_;
        end
        function g    = get.independentData(this)
            g = this.kernel_.independentData;
        end
        function this = set.independentData(this, s)
            this.kernel_.independentData = s;
        end
        function g    = get.dependentData(this)
            g = this.kernel_.dependentData;
        end
        function this = set.dependentData(this, s)
            this.kernel_.dependentData = s;
        end
        
        %%

        function ed   = estimateData(this)
            ed = this.kernel_.estimateData(this.modelParameters);
        end
        function Q    = objectiveFunc(this)
            Q = this.kernel_.objectiveFunc(this.modelParameters);
        end
        function ps   = modelParameters(varargin)
            ps = [];
        end
        function sps  = modelStdParameters(varargin)
            sps = [];
        end
        function        constructSyntheticKernel(~)
            error('mlbayesian:notImplemented', 'NullModel.constructSyntheticKernel');
        end
        function        constructKernelWithData(~)
            error('mlbayesian:notImplemented', 'NullModel.constructKernelWithData');
        end
        function ps   = solverParameters(this)
            switch (this.classOfSolver)
                case 'mlbayesian.LevenbergMarquardt'
                    ps = this.LMParameters;
                case 'mlbayesian.BretthorstMcmc'
                    ps = this.mcmcParameters;
                case 'mlnest.NestedSamplingMain'
                    ps = this.nestParameters;
                case 'mlstan.FlatHMC'
                    ps = this.hmcParameters;
                case 'mlstan.HierarchicalHMC'
                    ps = this.hierarchicalHmcParameters;
                otherwise
                    error('mlbayesian:unsupportedSwitchStrategy', ...
                        'PolynomialsModel.solverParameters.this.classOfSolver->%s', this.classOfSolver);
            end
        end	
        function this = updateModel(this, solvr)
            assert(isa(solvr, 'mlbayesian.ISolver'));
            switch (class(solvr))
                case 'mlbayesian.LevenbergMarquardt'
                    this = this.LMParameters2model(solvr);
                case 'mlbayesian.BretthorstMcmc'
                    this = this.mcmcParameters2model(solvr);
                case 'mlnest.NestedSamplingMain'
                    this = this.nestParameters2model(solvr);
                case 'mlstan.FlatHMC'
                    this = this.hmcParameters2model(solvr);
                case 'mlstan.HierarchicalHMC'
                    this = this.hierarchicalHmcParameters2model(solvr);
                otherwise
                    error('mlbayesian:unsupportedSwitchStrategy', ...
                        'PolynomialsModel.updateModel.solvr->%s', solvr);
            end
        end      

        function plot(this, varargin)
            figure;
            plot(this.independentData, this.dependentData, 'o',  ...
                 this.independentData, this.estimateData,  '-', varargin{:});
            legend('data', 'estimated');  
            title(class(this));
        end
        function writetable(varargin)
        end
		  
 		function this = NullModel(varargin)
 			%% NULLMODEL
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'independentData', [], @isnumeric);
            addParameter(ip, 'dependentData',   [], @isnumeric);
            addParameter(ip, 'solverClass',         @ischar);
            addParameter(ip, 'useSynthetic', false, @islogical);
            addParameter(ip, 'datedFilename', false, @islogical);
            parse(ip, varargin{:});  

            this.kernel_        = mlbayesian.NullKernel(ip.Results.independentData, ip.Results.dependentData);
            this.classOfSolver  = ip.Results.solverClass;
            this.useSynthetic   = ip.Results.useSynthetic;    
            this.datedFilename_ = ip.Results.datedFilename;
 		end
 	end 
    
    %% PROTECTED
    
    properties (Access = protected)
        kernel_
        datedFilename_
    end 
    
    methods (Access = protected)
        function this = setupKernel(this)
            if (this.useSynthetic || ...
                    all(isempty(this.independentData) || all(isempty(this.dependentData))))
                this = this.constructSyntheticKernel;
            else 
                this = this.constructKernelWithData;
            end  
        end
        function this = setupFilesystem(this)
            %% for mlio.AbstractIO            
            this.filepath_ = pwd;
            this.fileprefix_ = strrep(class(this), '.', '_');
            if (this.datedFilename_)
                this.fileprefix_ = [this.fileprefix_ '_' datestr(now, 30)];
            end
            this.filesuffix_ = '.mat';
        end
        function this = checkModel(this)
            assert(ischar(this.parameterIndexToName(length(this.modelParameters))), ... 
                'mismatched lengths of parameterIndexToName and modelParameters'); 
        end
        function ps   = LMParameters(this)
            ps = [];
        end
        function this = LMParameters2model(this, ~)
        end 
        function ps   = mcmcParameters(this)
            %% MCMCPARAMETERS must be in heap memory for speed
            %  @return struct containing:
            %  fixed      is logical, length := length(this.modelParameters)
            %  fixedValue is numeric, "
            %  min        is numeric, "; for prior distribution
            %  mean       is numeric, "
            %  max        is numeric, "; for prior distribution
            %  std        is numeric, "; for annealing
            %  nAnneal    =  20, number of loops per annealing temp
            %  nBeta      =  50, number of temperature steps
            %  nPop       =  50, number of population for annealing/burn-in and proposal/sampling
            %  nProposals = 100, number of proposals for importance sampling
            %  nSamples   is numeric, numel of independentData
            
            ps = struct( ...
                'fixed',      logical([]), ...
                'fixedValue', [], ...
                'min_',       [], ...
                'mean_',      [], ...
                'max_',       [], ... 
                'std_',       [], ...
                'nAnneal',    this.nAnneal, ...
                'nBeta',      this.nBeta, ...
                'nPop',       50, ...
                'nProposals', this.nProposals, ...
                'nSamples',   numel(this.independentData));
        end
        function this = mcmcParameters2model(this, solvr)
            assert(isa(solvr, 'mlbayesian.IMcmcSolver'));
            this.as  = ensureRowVector(solvr.kernel.bestFitParams);
            this.sas = ensureRowVector(solvr.kernel.stdParams);
        end
        function ps   = nestParameters(this)
            ps = [];
        end
        function this = nestParameters2model(this, ~)
        end 
        function ps   = hmcParameters(this)
            ps = [];
        end
        function this = hmcParameters2model(this, ~)
        end 
        function ps   = hierarchicalHmcParameters(this)
            ps = [];
        end
        function this = hierarchicalHmcParameters2model(this, ~)
        end 
    end  

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
    
 end

