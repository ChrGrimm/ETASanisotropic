function [ F, dF_D, dF_q ] = compute_infiniteSpatialIntegr( Catalog, ...
                                                             isCategory, ...
                                                             spatParamETAS )
% Spatial integral is approximated by 1, gradients by 0
% (if spatial kernel is not normalized to probability
% distribution, compute integral until event-specific spatial extent)
% Events outside of polygon are assumed to have spatial integral and gradiants equal to 0

    %% Initialize output vectors
    Catalog             = Catalog(isCategory,:);
    [ F, dF_D, dF_q ]   = initialize_spatIntegrOutputVectors( Catalog, 'default' );

    %% Input processing
    % Identify events in spatial target window
    isSpatTarget    = abs(Catalog.flag) >= 0.5;

    %% Compute spatial integral & gradiants for spatialtarget events
    % Since normalized, integral = 1, gradiants = 0
    F(isSpatTarget)     = 1;
    dF_D(isSpatTarget)  = 0;
    dF_q(isSpatTarget)  = 0;
    
end

%% Old code pieces
% % If not normalized, integrate till spatial extent
% f_inn_R             = f_inn(R, R.^2, magn, rupExtent, wi);
% f_inn_R_to1minusQ   = f_inn_R.^(1-q);
% F(isSpatTarget)     = 1 - f_inn_R_to1minusQ;
% dF_D(isSpatTarget)  = (1-q)/D * f_inn_R_to1minusQ .* (1-f_inn_R.^(-1)) ...
%                         .* isWiderThanMinWidth(isSpatTarget);
% dF_q(isSpatTarget)  = log(f_inn_R) .* f_inn_R_to1minusQ ...
%                         .* isWiderThanMinWidth(isSpatTarget);
