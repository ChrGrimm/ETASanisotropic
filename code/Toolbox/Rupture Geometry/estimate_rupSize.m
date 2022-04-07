function [ rupLength, rupWidth ] = estimate_rupSize( Mw, ...
                                                     tecZoneType, ...
                                                     faultingStyle, ...
                                                     spaceUnit )
                                                     %

    % tecZoneType should be vector - in case it is not, quantify
    if ~iscell(tecZoneType)
        tmp         = cell(size(Mw));
        tmp(:)      = {tecZoneType};
        tecZoneType = tmp;
    end
    % faultingStyle should be vector - in case it is not, quantify
    if ~iscell(faultingStyle)
        tmp         = cell(size(Mw));
        tmp(:)      = {faultingStyle};
        faultingStyle   = tmp;
    end
    
    % Initialization
    rupLength       = zeros(size(Mw));
    rupWidth        = zeros(size(Mw));
    
    % Marker for tectonic zone type
    isSubduction    = strcmp( tecZoneType, 'subduction' );
    isContinental   = strcmp( tecZoneType, 'continental' );
    
    % Marker for fault type
    isDummy         = strcmp( faultingStyle, 'dummy' );
    isReverse       = strcmp( faultingStyle, 'reverse' );
    isNormal        = strcmp( faultingStyle, 'normal' );
    isStrikeSlip    = strcmp( faultingStyle, 'strike-slip' );
    isAll           = strcmp( faultingStyle, 'all' );

    % Checks
    if ~ ( sum( isSubduction+isContinental ) == length(Mw) )
        error('Invalid tectonic zone type submitted!')
    end
    if ~ ( sum( isDummy+isReverse+isNormal+isStrikeSlip+isAll ) == length(Mw) )
        error('Invalid fault type submitted!')
    end
    
    %% Rupture length and width estimates
    % Subduction zone events: Blaser et al; Continental zone events: Wells & Coppersmith
    % dummy
    rupLength( isDummy )                    = 10^(-6) * ones(sum(isDummy),1);
    rupWidth( isDummy )                     = ones(sum(isDummy),1);
    % reverse
    rupLength( isReverse & isSubduction )   = 10.^( -2.37 + 0.57*Mw(isReverse & isSubduction) );
    rupWidth( isReverse & isSubduction )    = 10.^( -1.86 + 0.46*Mw(isReverse & isSubduction) );
    rupLength( isReverse & isContinental )  = 10.^( -2.42 + 0.58*Mw(isReverse & isContinental) );
    rupWidth( isReverse & isContinental )   = 10.^( -1.61 + 0.41*Mw(isReverse & isContinental) );
    % normal
    rupLength( isNormal & isSubduction )    = 10.^( -1.91 + 0.52*Mw(isNormal & isSubduction) );
    rupWidth( isNormal & isSubduction )     = 10.^( -1.20 + 0.36*Mw(isNormal & isSubduction) );
    rupLength( isNormal & isContinental )   = 10.^( -1.88 + 0.50*Mw(isNormal & isContinental) );
    rupWidth( isNormal & isContinental )    = 10.^( -1.14 + 0.35*Mw(isNormal & isContinental) );
    % strike-slip
    rupLength(isStrikeSlip & isSubduction)  = 10.^( -2.69 + 0.64*Mw(isStrikeSlip & isSubduction) );
    rupWidth(isStrikeSlip & isSubduction)   = 10.^( -1.12 + 0.33*Mw(isStrikeSlip & isSubduction) );
    rupLength(isStrikeSlip & isContinental) = 10.^( -2.57 + 0.62*Mw(isStrikeSlip & isContinental) );
    rupWidth(isStrikeSlip & isContinental)  = 10.^( -0.76 + 0.27*Mw(isStrikeSlip & isContinental) );
    % all
    rupLength( isAll & isSubduction )       = 10.^( -2.31 + 0.57*Mw(isAll & isSubduction) );
    rupWidth( isAll & isSubduction )        = 10.^( -1.56 + 0.41*Mw(isAll & isSubduction) );
    rupLength( isAll & isContinental )      = 10.^( -2.44 + 0.59*Mw(isAll & isContinental) );
    rupWidth( isAll & isContinental )       = 10.^( -1.01 + 0.32*Mw(isAll & isContinental) );
    
    %% Optional: Conversion to degree units
    if strcmp( spaceUnit, 'degree' )
        rupLength   = rupLength/111;
        rupWidth    = rupWidth/111;
    end
    
%     Kagan and Jackson (1999) 
% 
%         % formula (1) for converting GEM moment magnitude formate to scalar
%         % seismic moment M0
%         M0 = 10.^(3/2*(Mw + 6.0));
% 
%         % According to Kagan-Jackson (1999) below formula (2)
%         W = 25000; % in meters (m), rupture width
%         dS = 30*10^5; % in N/m^2 (=30 bar), earthquake stress drop
% 
%         % formula (3) for estimating rupture length (in km)
%         rupLength = sqrt(M0/W/dS)/1000;


end

