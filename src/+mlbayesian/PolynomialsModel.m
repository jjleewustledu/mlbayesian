classdef PolynomialsModel < mlanalysis.NullModel
	%% POLYNOMIALSMODEL supports p(t) := \Sigma_{i = 0}^3 a_i t^i

	%  $Revision$
 	%  was created 24-Dec-2017 17:30:19 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    
    properties (Dependent)        
        asMin
        asMax
    end
    
	properties
        t0 = 0;
        dt = 1;
        tfinal = 10;
        as = [1; .1; .01; .001] % expected polynomial coefficients
        sas = [1; .1; .01; .001] % std used for annealing        
        fixed = false(4,1)
        fixedValue = nan(4,1)
    end
    
    methods (Static)
        function name = parameterIndexToName(idx)
            names = {'a1' 'a2' 'a3' 'a4'}; % per Matlab's Fortran indexing
            name = names{idx};
        end
    end

	methods 
        
        %% GET, SET
        
        function g = get.asMin(this)
            if (~isempty(this.asMin_))
                g = this.asMin_;
                return
            end
            g = this.as - this.M*this.sas;
            g = ensureColVector(g);
        end
        function this = set.asMin(this, s)
            assert(isnumeric(s));
            assert(numel(this.as) == numel(s));
            this.asMin_ = s;
        end
        function g = get.asMax(this)
            if (~isempty(this.asMax_))
                g = this.asMax_;
                return
            end
            g = this.as + this.M*this.sas;
            g = ensureColVector(g);
        end
        function this = set.asMax(this, s)
            assert(isnumeric(s));
            assert(numel(this.as) == numel(s));
            this.asMax_ = s;
        end
                
        %%
		  
        function ps   = modelParameters(this)
            ps = ensureColVector(this.as);
        end
        function sps  = modelStdParameters(this)
            sps = ensureColVector(this.sas);
        end
        function this = doConstructGenerative(this)
            idata = this.t0:this.dt:this.tfinal;  
            ddata = mlbayesian.PolynomialsKernel.polynomial(this.as, idata);
            this.kernel_ = mlbayesian.PolynomialsKernel(idata, ddata);
        end
        
 		function this = PolynomialsModel(varargin)
 			%% POLYNOMIALSMODEL
            %  @param independentData is numeric, defaults to this.t0:this.dt:this.tfinal.
            %  @param dependentData   is numeric, defaults to generative model.
            %  @param constructGenerative is logical; forces dependentData := generative model if true.
            %  @returns mlanalysis.IModel solvable by mlanalysis.ISolver implementations.
            
            this = this@mlanalysis.NullModel(varargin{:});            
            if (this.constructGenerative || ...
                    all(isempty(this.independentData) || all(isempty(this.dependentData))))
                this = this.doConstructGenerative;
            else                          
                this.kernel_ = mlbayesian.PolynomialsKernel( ...
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
    
    properties (Access = protected)
        asMin_
        asMax_
    end
    
    methods (Access = protected)
        function ps = mcmcParameters(this)
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
                'fixed',      ensureColVector(this.fixed), ...
                'fixedValue', ensureColVector(this.fixedValue), ...
                'min_',       this.asMin, ...
                'mean_',      this.modelParameters, ...
                'max_',       this.asMax, ... 
                'std_',       this.modelStdParameters, ...
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
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

