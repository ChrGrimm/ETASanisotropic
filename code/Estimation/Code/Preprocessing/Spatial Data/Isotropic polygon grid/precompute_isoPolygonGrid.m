function [ SpatData, ...
           typeKernel ] = precompute_isoPolygonGrid( SpatData, ...
                                                     polygon, ...
                                                     Catalog, ...
                                                     spaceModel, ...
                                                     spaceUnit, ...
                                                     makeChecks, ...
                                                     makePlots )
% 
% This function precomputes the spatial integral over the ETAS target space window of the isotropic 
% Gaussian kernel density used to redistribute the location of background events.
% If modelDesign is chosen to be isotropic, the related information for computing spatial integrals 
% of isotropic spatial triggering distributions in later ETAS iteration steps are stored in cell
% arrays.
%
% Call: [ backgrIntegral, ...
%           r_squ, ...
%           r0_squ, ...
%           segmWeight, ...
%           isLower, ...
%           isUpper ] = precompute_isoSpatDistr( polygon, Catalog, bandwidth, modelDesign, spaceUnit )
%
% INPUTS:
%
% - polygon:           numeric array (Px2)     (centralized) polygon coordinates (x | y)
% - Catalog:        table                   historic event set. Relevant columns in this function:
%                                           - x: event coordinate x
%                                           - y: event coordinate y
%                                           - flag: determines type of event: >= 0 for target time
% - bandwidth:      numeric vector (Bx1)    bandwidths value for all events with flag >= 0
% - modelDesign:    string                  indicates the chosen model design (e.g. isotropic, ...)
% - spaceUnit:      string                  indicates the chosen space unit (e.g. km, degrees, ...)
%
% OUTPUTS:
%
% - backgrIntegral: numeric vector (Bx1)    background integral for all events with flag >= 0
% - r_squ:          cell vector (numerics)  squared event-to-polygon distances for polygon grid
%                                           points
% - r0_squ:         cell vector (numerics)  squared event-to-polygon distances for intermediate
%                                           polygon grid points
% - segmWeight:     cell vector (numerics)  weight of the respective radial segment to entire circle
% - isLower:        cell vector (booleans)  booleans indicating if r_squ values are lower boundary
% - isUpper:        cell vector (booleans)  booleans indicating if r_squ values are upper boundary

    
    idxIso      = find(strcmp(Catalog.typeKernel, 'iso') | Catalog.flag >=0);
    typeKernel  = Catalog.typeKernel;
    
    isTargetT   = Catalog.flag >= 0 & ~Catalog.isDupl;
    
    %% Initializations
    % cell arrays
    SpatData = initialize_tables4precomputations( SpatData, length(Catalog.x) );
    
    % start and end points of polygon edges
    polyStart   = polygon(1:end-1,:);
    polyEnd     = polygon(2:end,:);
    
    % background integrals
    SpatData.backgrIntegral  = zeros(sum(isTargetT),1);
    
    % target time event index
    jTargetEvent    = 1;
    
    %% Computation of radial segments and background integrals
    % Loop over all events
    for iEvent = idxIso'
        
        [ r, ...
          r0, ...
          r_squ, ...
          r0_squ, ...
          segmFactor, ...
          isStartP, ...
          isEndP ] = compute_isotropicGrid( polyStart, ...
                                             polyEnd, ...
                                             iEvent, ...
                                             Catalog.x(iEvent), ...
                                             Catalog.y(iEvent), ...
                                             Catalog.flag(iEvent), ...
                                             spaceUnit, ...
                                             'isotropic', ...
                                             makeChecks, ...
                                             makePlots );

%             % Check if sum of segment weights and factors equal desired value
%             check_isoSpatDistr( IsoGrid, Catalog.flag(iEvent), Catalog, iEvent )

        %% Background integral
        % Compute background integral for target time events only            
        if isTargetT(iEvent)

            if strcmp( spaceModel, 'none')
                SpatData.backgrIntegral(jTargetEvent)   = Catalog.flag(iEvent) > 0;
            else
                SpatData.backgrIntegral(jTargetEvent) = compute_backgrIntegral( SpatData.bandwidth(jTargetEvent), ....
                                                                                r_squ, ...
                                                                                r0_squ, ...
                                                                                segmFactor, ...
                                                                                isStartP, ...
                                                                                isEndP );
            end
            jTargetEvent                            = jTargetEvent + 1;
            
        end

        %% Store data in tables if model is run in iso-mode
        if strcmp(Catalog.typeKernel{iEvent}, 'iso')

            [ r, r_squ, r0, r0_squ, segmFactor, isStartP, isEndP ] ...
                = apply_spatialExtent2SpatData( Catalog.spatRestr(iEvent), ...
                                                r, r0, segmFactor, isStartP, isEndP );

            % Store event-to-polygon distances, segment factors and start-end booleans
            SpatData.Iso_r{iEvent}              = r;
            SpatData.Iso_r0{iEvent}             = r0;
            SpatData.Iso_r_squ{iEvent}          = r_squ;
            SpatData.Iso_r0_squ{iEvent}         = r0_squ;
            SpatData.Iso_segmFactor{iEvent}     = segmFactor;
            SpatData.Iso_isStartP{iEvent}       = isStartP;
            SpatData.Iso_isEndP{iEvent}         = isEndP;
            
            % Change type of kernel if event is completely inside or completely outside of polygon
            if length(r)==1
                typeKernel{iEvent}      = 'full';
            end
            
        end
            
%         end       
        
    end

end
