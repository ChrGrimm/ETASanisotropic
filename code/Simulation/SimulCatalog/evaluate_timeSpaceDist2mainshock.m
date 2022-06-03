function [ tDaysDiff, ...
           distKm ] = evaluate_timeSpaceDist2mainshock( SynthCatalog )
       
    tDaysDiff   = -1*ones(size(SynthCatalog,1),1);
    distKm      = -1*ones(size(SynthCatalog,1),1);
    idxTriggered= find( SynthCatalog.flag ~= 1 );
    
    idxMainsh   = SynthCatalog.idxMainsh(idxTriggered);
    
    tDaysDiff(idxTriggered) = SynthCatalog.t(idxTriggered) - SynthCatalog.t(idxMainsh);
    
    [ rupStart, rupEnd ] = compute_rupCoords( SynthCatalog.x(idxMainsh), ...
                                              SynthCatalog.y(idxMainsh), ...
                                              SynthCatalog.rupExtent(idxMainsh), ...
                                              SynthCatalog.strike(idxMainsh), ...
                                              SynthCatalog.epiPos(idxMainsh) );
    
    for i=1:length(idxTriggered)
        distKm(idxTriggered(i)) = 111 * get_dist2line( rupStart(i,:), rupEnd(i,:), ...
                                                       [SynthCatalog.x(idxTriggered(i)), SynthCatalog.y(idxTriggered(i))]', ...
                                                       'degree', 0 );   
    end
       
end