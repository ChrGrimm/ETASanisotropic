function [xx_new, yy_new] = sample_spaceBackground( nEvents, BackgrGrid )

    %% Compute cumulative sum of all rates to sample from them
    idxGridP    = randsample( 1:length(BackgrGrid.gridX), nEvents, true, BackgrGrid.distr );    
    xx_new      = BackgrGrid.gridX(idxGridP) + BackgrGrid.deltaX * (rand(nEvents,1)-0.5);
    yy_new      = BackgrGrid.gridY(idxGridP) + BackgrGrid.deltaY * (rand(nEvents,1)-0.5);

end

