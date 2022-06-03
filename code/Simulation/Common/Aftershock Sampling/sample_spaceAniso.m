function [ xOffspring, yOffspring ] = sample_spaceAniso( xTrigger, yTrigger, rupLTrigger, strikeTrigger, epiPosTrigger, r )

    xOffspring  = zeros(size(xTrigger));
    yOffspring  = zeros(size(yTrigger));
    u           = rand( length(xTrigger), 1 );

    [ rupStart, ...
      rupEnd ] = compute_rupCoords( xTrigger, yTrigger, rupLTrigger, strikeTrigger, epiPosTrigger );

    contourPoint = (2*pi*r + 2*rupLTrigger) .* u;

    isLeftCircle    = contourPoint < pi*r;
    isRightCircle   = contourPoint >= pi*r & contourPoint < 2*pi*r;
    isAboveLine     = contourPoint >= 2*pi*r & contourPoint < 2*pi*r+rupLTrigger;
    isBelowLine     = contourPoint >= 2*pi*r+rupLTrigger;

    %% isLeftCircle
    alphaStrike                 = 90 - strikeTrigger(isLeftCircle);
    alphaAftershock             = deg2rad( 90+alphaStrike ) + contourPoint( isLeftCircle ) ./ r(isLeftCircle);
    xOffspring( isLeftCircle )  = rupStart(isLeftCircle,1) + r(isLeftCircle) .* cos(alphaAftershock);
    yOffspring( isLeftCircle )  = rupStart(isLeftCircle,2) + r(isLeftCircle) .* sin(alphaAftershock);

    %% isRightCircle
    alphaStrike                 = 90 - strikeTrigger(isRightCircle);
    alphaAftershock             = deg2rad( 90+alphaStrike ) - ( contourPoint( isRightCircle ) ./ r(isRightCircle) - pi );
    xOffspring( isRightCircle ) = rupEnd(isRightCircle,1) + r(isRightCircle) .* cos(alphaAftershock);
    yOffspring( isRightCircle ) = rupEnd(isRightCircle,2) + r(isRightCircle) .* sin(alphaAftershock);

    %% isAboveLine        
    m                       = rupEnd(isAboveLine,:) - rupStart(isAboveLine,:);              % directional vector
    proj                    = ( contourPoint(isAboveLine)-2*pi*r(isAboveLine) ) ./ rupLTrigger(isAboveLine);
    if isempty(proj)  
        proj = ones(0,2);  
    end
    refPointOnLine          = rupStart(isAboveLine,:) + proj .* m;
    alphaStrike             = 90 - strikeTrigger(isAboveLine);
    alphaAftershock         = deg2rad( 90+alphaStrike );
    xOffspring(isAboveLine) = refPointOnLine(:,1) + r(isAboveLine) .* cos(alphaAftershock);
    yOffspring(isAboveLine) = refPointOnLine(:,2) + r(isAboveLine) .* sin(alphaAftershock);

   %% isBelowLine        
    m                       = rupEnd(isBelowLine,:) - rupStart(isBelowLine,:);              % directional vector
    proj                    = ( contourPoint(isBelowLine) - ( 2*pi*r(isBelowLine)+rupLTrigger(isBelowLine) ) ) ./ rupLTrigger(isBelowLine);
    if isempty(proj)
        proj = ones(0,2);  
    end
    refPointOnLine          = rupStart(isBelowLine,:) + proj .* m;
    alphaStrike             = 90 - strikeTrigger(isBelowLine);
    alphaAftershock         = deg2rad( alphaStrike-90 );
    xOffspring(isBelowLine) = refPointOnLine(:,1) + r(isBelowLine) .* cos(alphaAftershock);
    yOffspring(isBelowLine) = refPointOnLine(:,2) + r(isBelowLine) .* sin(alphaAftershock);
    
end