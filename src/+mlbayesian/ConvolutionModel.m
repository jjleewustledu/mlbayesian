classdef ConvolutionModel < mlanalysis.NullModel
	%% CONVOLUTIONMODEL supports p(t) := \Sigma_{i = 0}^3 a_i t^i

	%  $Revision$
 	%  was created 24-Dec-2017 17:30:19 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    
	properties
        t0 = 0
        dt = 1
        tfinal = 10
        A = 1
        T = 1
        U = 1
        sA = 1
        sT = 1
        sU = 1
        
        fixed = false(3,1)
        fixedValue = nan(3,1)
    end
    
    methods (Static)
        function name = parameterIndexToName(idx)
            names = {'A' 'T' 'U'}; % per Matlab's Fortran indexing
            name = names{idx};
        end
    end

	methods 
		  
        function ps   = modelParameters(this)
            ps = [this.A; this.T; this.U];
        end
        function sps  = modelStdParameters(this)
            sps = [this.sA; this.sT; this.sU];
        end
        function this = doConstructGenerative(this)
            idata = this.t0:this.dt:this.tfinal; 
            ddata = mlbayesian.ConvolutionKernel.convolution(this.modelParameters, idata);
            this.kernel_ = mlbayesian.ConvolutionKernel(idata, ddata);
        end
        
 		function this = ConvolutionModel(varargin)
 			%% CONVOLUTIONMODEL
            %  @param independentData is numeric, defaults to this.t0:this.dt:this.tfinal.
            %  @param dependentData   is numeric, defaults to generative model.
            %  @param constructGenerative is logical; forces dependentData := generative model if true.
            %  @returns mlanalysis.IModel solvable by mlanalysis.ISolver implementations.
            
            this = this@mlanalysis.NullModel(varargin{:});
            if (this.constructGenerative || ...
                    all(isempty(this.independentData) || all(isempty(this.dependentData))))
                this = this.doConstructGenerative;
            else 
                this.kernel_ = mlbayesian.ConvolutionKernel( ...
                    this.independentData, this.dependentData);
            end            
            
            assert(ischar(this.parameterIndexToName(length(this.modelParameters))), ... 
                'mismatched lengths of parameterIndexToName and modelParameters');
            
            %% for mlio.AbstractIO
            
            this.filepath_ = pwd;
            this.fileprefix_ = strrep(class(this), '.', '_');
            if (this.datedFilename_)
                this.fileprefix_ = [this.fileprefix_ '_' datestr(now, 30)];
            end
            this.filesuffix_ = '.mat';
 		end
    end 
    
    %% PROTECTED  
    
    methods (Access = protected)
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
            
            ps = ensureColVector(this.modelParameters);
            sps = ensureColVector(this.modelStdParameters);
            ps = struct( ...
                'fixed',      ensureColVector(this.fixed), ...
                'fixedValue', ensureColVector(this.fixedValue), ...
                'min_',       ps - this.M*sps, ...
                'mean_',      ps, ...
                'max_',       ps + this.M*sps, ... 
                'std_',       sps, ...
                'nAnneal',    this.nAnneal, ...
                'nBeta',      this.nBeta, ...
                'nPop',       50, ...
                'nProposals', this.nProposals, ...
                'nSamples',   numel(this.independentData));
        end
        function this = mcmcParameters2model(this, solvr)
            assert(isa(solvr, 'mlbayesian.IMcmcSolver'));
            bfp = ensureRowVector(solvr.kernel.bestFitParams);
            sp  = ensureRowVector(solvr.kernel.stdParams);
            this.A = bfp(1);
            this.T = bfp(2);
            this.U = bfp(3);
            this.sA = sp(1);
            this.sT = sp(2);
            this.sU = sp(3);
        end
    end    

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

