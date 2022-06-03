function [ xOffspring, yOffspring ] = sample_spaceIso( xTrigger, yTrigger, r )

    % Sample angle [0, 360]
    phi = 2*pi*rand(length(xTrigger),1);

    % Evaluate movements dX and dY
    xOffspring = xTrigger + r .* cos(phi); % /cos(cntLat*pi/180)
    yOffspring = yTrigger + r .* sin(phi);
    
end