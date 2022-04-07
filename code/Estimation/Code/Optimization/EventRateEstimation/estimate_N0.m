function [ N0, ...
           dN0 ] = estimate_N0( Catalog, ...
                                   Inputs, ...
                                   TempData, ...
                                   nBackgrPerDay, ...
                                   sqrtParamETAS, ...
                                   ModelFuncs, ...
                                   F, dF_D, dF_gamma, dF_q, ...
                                   strMode )
% This function estimates the "true" temporal event rate N0 at a time grid covering the target time
% period. In the first part, event rate contributions from the most recent events are computed.
% Then, event rate contributions from other past events are approximated linearly between two
% consecutive time grid points.

    %% Extract information
    % Extract current parameter values (they are in square root format  for numerical reasons)
    [mu, A, ~, c, p, ~, ~, ~, Tb] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Extract event-specific column information from table Catalog
    t           = Catalog.t;
    mag         = Catalog.mag;
    evWeight    = Catalog.evWeight;
    isTargetTime= Catalog.flag >= 0;
    tempRestr   = Catalog.tempRestr;
    
    % Determine whether gradients should be computed
    computeGradients    = strcmp(strMode, 'with gradient');
    
    % Compute constant, event-specific factors
    Mc = Inputs.TargetSettings.Mc;
    productivity        = ModelFuncs.fProductivity( mag, Mc );
    factor              = productivity .* evWeight .* F;
    
    % Extract time grid data
    idxTimeGrid             = TempData.idxTimeGrid;
    idxRecentEvents         = TempData.idxRecentEvents;
    timeSpan4recentEvents   = TempData.timeSpan4recentEvents;
    tFixGrid                = TempData.tFixGrid;
    thresh4largeMag         = TempData.thresh4largeMag;
    
    fixedParamETAS          = Inputs.ModelSettings.fixedParamETAS;
    spaceModel              = Inputs.ModelSettings.spaceModel;
    Mc                      = Inputs.TargetSettings.Mc;
    
    %% Compute event rate contributions from recent events
    % Compute Omori-Utsu law term for time differences between grid points and events
    fTempK_inn_ij            = TempData.timeDiff + c;
    fTempK_inn_ij_toMinusP   = fTempK_inn_ij.^(-p);
    
    % Compute trigger rate contributions as a product of event-specific factor and Omori-Utsu term
    N0_trig_ij          = factor(idxRecentEvents) .* fTempK_inn_ij_toMinusP;
    % Sum up event contributions that belong to the same time grid point
    N0_trig_recent      = accumarray( idxTimeGrid, N0_trig_ij )';
    
    % Compute gradients of event rate contributions from recent events
    if computeGradients
        
        dN0_recent = zeros( length(sqrtParamETAS), length(N0_trig_recent) );
        
        if ~fixedParamETAS(2)
            dN0_recent(2,:) = Tb * N0_trig_recent/A;
        end
        if ~fixedParamETAS(3)
            dN0_recent(3,:) = Tb * accumarray( idxTimeGrid, N0_trig_ij .* ( mag(idxRecentEvents)-Mc ) );
        end
        if ~fixedParamETAS(4)
            tmp_c           = -p ./ fTempK_inn_ij;
            dN0_recent(4,:) = Tb * accumarray( idxTimeGrid, N0_trig_ij .* tmp_c );
        end
        if ~fixedParamETAS(5)
            tmp_p           = -log(fTempK_inn_ij);
            dN0_recent(5,:) = Tb * accumarray( idxTimeGrid, N0_trig_ij .* tmp_p );
        end
        if ~strcmp(spaceModel, 'none')
            
            % Extract specific TempData for spatial gradients
            idxTimeGrid_spatGrad            = TempData.idxTimeGrid_spatGrad;
            idxRecentEvents_spatGrad        = TempData.idxRecentEvents_spatGrad;
            fTempK_inn_ij_toMinusP_spatGrad = fTempK_inn_ij_toMinusP( ~TempData.isZeroSpatialGradient );
            
            factor_dD       = productivity .* evWeight .* dF_D;
            factor_dGamma   = productivity .* evWeight .* dF_gamma;
            factor_dQ       = productivity .* evWeight .* dF_q;
            if ~isempty(idxTimeGrid_spatGrad)
                idxCol = 1:max(idxTimeGrid_spatGrad);
                if ~fixedParamETAS(6)
                    dN0_recent(6,idxCol) = Tb * accumarray( idxTimeGrid_spatGrad, factor_dD(idxRecentEvents_spatGrad) .* fTempK_inn_ij_toMinusP_spatGrad );
                end
                if ~fixedParamETAS(7)
                    dN0_recent(7,idxCol) = Tb * accumarray( idxTimeGrid_spatGrad, factor_dGamma(idxRecentEvents_spatGrad) .* fTempK_inn_ij_toMinusP_spatGrad );
                end
                if ~fixedParamETAS(8)
                    dN0_recent(8,idxCol) = Tb * accumarray( idxTimeGrid_spatGrad, factor_dQ(idxRecentEvents_spatGrad) .* fTempK_inn_ij_toMinusP_spatGrad );
                end
            end    
        end
    end
    
%     % Add time-homogeneous background rate
%     N0_recent          = Tb * (N0_trig_recent + mu * nBackgrPerDay);
%     dN0_recent         = Tb * dN0_recent;
%     dN0_recent(9,:)    = N0_trig_recent/Tb;
%     dN0_recent         = dN0_recent * 2 .* sqrtParamETAS;
    
    %% Approximate event rate contributions from other events linearly
%     % Set fix dates (event dates, event dates - epsTime, start and end time)
%     epsTime         = 10^-10;
%     eventTimes      = [0; t(isTargetTime); timeWindow_end]';
%     tFixGrid        = sort([0; t(isTargetTime); t(isTargetTime)-epsTime; timeWindow_end]');
    
    % Initialize matrices to save results
    nTimeGrid           = length( tFixGrid );
    N0_trig_other       = zeros( 1, nTimeGrid );
    dN0_other           = zeros( length(sqrtParamETAS), nTimeGrid );
    isZeroSpatGradient  = strcmp( Catalog.typeKernel, 'full' );
    
    % Compute gradients of event rate contributions from other events
    for i = 1:nTimeGrid/2
        
        % indices of current grid points
        idxTimePoints   = [2*i-1, 2*i];
        tStart          = tFixGrid(2*i-1);

        % Evaluate which events are contributing, but are not recent
        isOtherEvents       = t < tStart - timeSpan4recentEvents(i) & ...
                              t >= tStart - tempRestr & ...
                              mag <= thresh4largeMag;
        if ~any(isOtherEvents); continue; end
        
        % Compute Omori-Utsu law term at start and end point of current time interval
        fTempK_inn_ij            = tFixGrid(idxTimePoints) - t(isOtherEvents) + c;
        fTempK_inn_ij_toMinusP   = fTempK_inn_ij.^(-p); 
                
        % Compute trigger rate contributions as a product of event-specific factor and Omori-Utsu term
        N0_trig_ij          = factor(isOtherEvents) .*  fTempK_inn_ij_toMinusP;
        
        % Plausibility check
        if any(~isreal(N0_trig_ij)) 
            error('Imaginary result!')
        end
        
        % Sum up event rate contributions
        N0_trig_other(idxTimePoints) = sum( N0_trig_ij, 1 );
                           
        %% Computation of temporal and spatial event rate
        if computeGradients  
            
            % Gradiants by productivity parameters
            if ~fixedParamETAS(2)
                dN0_other(2,idxTimePoints)    = Tb * N0_trig_other(idxTimePoints)/A;
            end
            if ~fixedParamETAS(3)
                dN0_other(3,idxTimePoints)    = Tb * sum( N0_trig_ij .* (mag(isOtherEvents)-Mc), 1 );
            end
            if ~fixedParamETAS(4)
                tmp_c                       = -p ./ fTempK_inn_ij;
                dN0_other(4,idxTimePoints)    = Tb * sum( N0_trig_ij .* tmp_c, 1 );
            end
            if ~fixedParamETAS(5)
                tmp_p                       = -log(fTempK_inn_ij);
                dN0_other(5,idxTimePoints)    = Tb * sum( N0_trig_ij .* tmp_p, 1 );
            end
            if ~strcmp(spaceModel, 'none')
                isOtherEvents = isOtherEvents & ~isZeroSpatGradient;
                isNoZeroGrad  = ~isZeroSpatGradient( isOtherEvents );
                if ~fixedParamETAS(6)
                    dN0_other(6,idxTimePoints)    = Tb * sum( factor_dD(isOtherEvents) .* fTempK_inn_ij_toMinusP(isNoZeroGrad,:), 1 );
                end
                if ~fixedParamETAS(7)
                    dN0_other(7,idxTimePoints)    = Tb * sum( factor_dGamma(isOtherEvents) .* fTempK_inn_ij_toMinusP(isNoZeroGrad,:), 1 );
                end
                if ~fixedParamETAS(8)
                    dN0_other(8,idxTimePoints)    = Tb * sum( factor_dQ(isOtherEvents) .* fTempK_inn_ij_toMinusP(isNoZeroGrad,:), 1 );
                end
            end
        end
        
    end
    
    %% Interpolation
    x               = tFixGrid;
    newX            = TempData.timeGrid';
    N0_trig_other   = interp1(x, N0_trig_other, newX);
    dN0_other       = interp1(x', dN0_other', newX)';
    
    % Add time-homogeneous background rate
    N0          = Tb * (  mu * nBackgrPerDay + N0_trig_recent + N0_trig_other );
    
    % Gradients
    if computeGradients
        dN0     = dN0_recent + dN0_other;
        if ~fixedParamETAS(1)
            if strcmp( spaceModel, 'none' )
                dN0(1,:)            = Tb;
            else
                dN0(1,:)            = Tb * nBackgrPerDay;
            end
        end
        dN0         = dN0 * 2 .* sqrtParamETAS;
        if ~fixedParamETAS(9)
            dN0(9,:)    = N0/Tb;
        end
    end
    
end
