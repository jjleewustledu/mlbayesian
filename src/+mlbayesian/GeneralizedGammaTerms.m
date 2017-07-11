classdef GeneralizedGammaTerms 
	%% GENERALIZEDGAMMATERMS provides static function for gamma distributions and generalization optimized for flops.
    %  There are no checks on input validity.  
    %  https://en.wikipedia.org/wiki/Gamma_distribution#Characterization_using_shape_.CE.B1_and_rate_.CE.B2
    %  https://en.wikipedia.org/wiki/Generalized_gamma_distribution
    %  N.B.  \rho(t; a,b,p) = \frac{p b^a}{\Gamma(a/p)} t^{a-1} e^{-(b t)^p} \text{ with } 1/b > 0, a > 0, p > 0, t > 0.

	%  $Revision$
 	%  was created 27-Jun-2017 22:38:01 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlbayesian/src/+mlbayesian.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	methods (Static)
        function rho = gammaStretchSeriesSteady(a1, b1, p1, t01, a2, b2, p2, t02, weight1, S, k, t0, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = steadyState(S, k, t0, t) + ...
                weight1*gammaStretchC(a1, b1, p1, t01, t) + (1 - weight1)*gammaStretchC(a2, b2, p2, t02, t);
            rho = abs(rho);
        end      
        function rho = gammaStretchSeries(a1, b1, p1, t01, a2, b2, p2, t02, weight1, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = weight1*gammaStretchC(a1, b1, p1, t01, t) + (1 - weight1)*gammaStretchC(a2, b2, p2, t02, t);
            rho = abs(rho);
        end    
        function rho = gammaStretchPair(a1, b1, p1, t01, t02, weight1, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = weight1*gammaStretchC(a1, b1, p1, t01, t) + (1 - weight1)*gammaStretchC(a1, b1, p1, t02, t);
            rho = abs(rho);
        end        
        function rho = gammaSeriesSteady(a1, b1, t01, a2, b2, t02, weight1, S, k, t0, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = steadyState(S, k, t0, t) + ...
                weight1*gammaTermR(a1, b1, t01, t) + (1 - weight1)*gammaTermR(a2, b2, t02, t);
        end        
        function rho = gammaSeries(a1, b1, t01, a2, b2, t02, weight1, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = weight1*gammaTermR(a1, b1, t01, t) + (1 - weight1)*gammaTermR(a2, b2, t02, t);
        end
        function rho = gammaPair(a1, b1, t01, t02, weight1, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = weight1*gammaTermR(a1, b1, t01, t) + (1 - weight1)*gammaTermR(a1, b1, t02, t);
        end
        function rho = gammaStretchSteady(a1, b1, p1, t01, S, k, t0, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = steadyState(S, k, t0, t) + gammaStretchC(a1, b1, p1, t01, t);
            rho = abs(rho);
        end    
        function rho = gammaStretch(a, b, p, t0, t)
            %% GAMMASTRETCH
            %  @returns \rho(t; a,b,p) = \frac{p b^a}{\Gamma(a/p)} t^{a-1} e^{-(b t)^p} \text{ with } 1/b > 0, a > 0, p > 0, t > 0.
            
            rho = abs(mlbayesian.GeneralizedGammaTerms.gammaStretchC(a, b, p, t0, t));
        end
        function rho = gammaSteady(a1, b1, t01, S, k, t0, t)
            import mlbayesian.GeneralizedGammaTerms.*;
            rho = steadyState(S, k, t0, t) + gammaTermR(a1, b1, t01, t);
        end
        function rho = gammaTerm(a, b, t0, t)
            %% GAMMATERM
            %  @returns \rho(t; a,b,p) = \frac{b^a}{\Gamma(a)} t^{a-1} e^{-b t} \text{ with } 1/b > 0, a > 0, t > 0.
            
            rho = mlbayesian.GeneralizedGammaTerms.gammaTermR(a, b, t0, t);
        end
        function rho = steadyState(S, k, t0, t)
            %% STEADYSTATE
            %  @returns \rho(t; S,k) =S(1 - e^{-k t}) \text{Heaviside}(t, t_0)
            
            if (t(1) >= t0) % saves extra flops from slide()
                tau = t - t0;
                rho = S*(1 - exp(-k*tau)) .* mlbayesian.AbstractBayesianStrategy.Heaviside(tau, t0);
            else 
                tau = t - t(1);
                rho = S*(1 - exp(-k*tau));
                rho = mlbayesian.AbstractBayesianStrategy.slide(rho, t, t0 - t(1));
            end
        end
    end 
    
    methods (Static, Access = protected)
        
        %% designed for flops efficiency
        
        function rho = gammaStretchC(a, b, p, t0, t)
            if (t(1) >= t0) % saves extra flops from slide(), \int dt rho(t) < 1
                tau = t - t0;
                rho = tau.^(a-1) .* exp(-(b*tau).^p);
            else
                tau = t - t(1);
                rho = tau.^(a-1) .* exp(-(b*tau).^p);
                rho = mlbayesian.AbstractBayesianStrategy.slide(rho, t, t0 - t(1));
            end
            rho = rho/sum(abs(rho));
            %rho = rho*p*b^a/gamma(a/p);
        end
        function rho = gammaTermR(a, b, t0, t)
            if (t(1) >= t0) % saves extra flops from slide(), \int dt rho(t) < 1
                tau = t - t0;
                rho = tau.^(a-1) .* exp(-b*tau);
            else
                tau = t - t(1);
                rho = tau.^(a-1) .* exp(-b*tau);
                rho = mlbayesian.AbstractBayesianStrategy.slide(rho, t, t0 - t(1));
            end
            rho = rho/sum(rho);
            %rho = rho*b^a/gamma(a);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

