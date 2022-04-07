function BranchingRatio = compute_branchingRatio( Catalog, ...
                                                    isTarget, ...
                                                    beta, ...
                                                    TargetSettings, ...
                                                    ModelFuncs, ...
                                                    version )
%
% This function computes the branching ratio according to the estimated ETAS model parameters.
%
% Call: branchRatio = compute_branchingRatio( paramETAS, beta )
%
% INPUTS:
%
% - beta:           numeric vector (2x1)    vector of Gutenberg-Richter parameter beta estimates:
%                                           (log-likelihood, root mean square)
%
% OUTPUTS:
%
% - branchRatio     numeric     branching ratio of estimated ETAS model based on:
%                                           log-likelihood
        
    %% Definition 1: Integral over aftershock productivity * Gutenberg-Richter
    Mc                  = TargetSettings.Mc;
    Mmax                = TargetSettings.Mmax;
    fMagnDistr          = @(mag) beta*exp(-beta*(mag-Mc)) / (1-exp(-beta*(Mmax-Mc)));
    f_branching         = @(mag) ModelFuncs.fProductivity_scaled(mag,Mc) .* fMagnDistr(mag);
    BranchingRatio(1)   = integral(f_branching, Mc, Mmax);
    
%     if strcmp(version, 'estimation')
%         %% Definition 2: Average probability of being triggered event
%         BranchingRatio(2)   = mean(1-Catalog.backgrProb(isTarget));
% 
%     end

end
