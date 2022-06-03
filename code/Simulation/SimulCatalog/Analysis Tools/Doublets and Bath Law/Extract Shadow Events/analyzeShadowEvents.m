ShadowEvents = table;
iRow = 0;
warning off
for i=1:length(idxEvents)
    if ~isempty(shadowingEvents{i})
        %
        iRow = iRow + 1;
                %
        idxRefEvents                        = shadowingEvents{i};
        idxLatest                           = max(idxRefEvents);
        [~,idxStrongest]                    = max( Catalog.mag(idxRefEvents) );
        idxStrongest                        = idxRefEvents( idxStrongest );
        %
        ShadowEvents.id(iRow)               = Catalog.id( idxEvents(i) );
        ShadowEvents.idsBefore{iRow}        = Catalog.id( shadowingEvents{i} );
        ShadowEvents.idLatest(iRow)         = Catalog.id( idxLatest );
        ShadowEvents.idStrongest(iRow)      = Catalog.id( idxStrongest );
        ShadowEvents.date(iRow)             = Catalog.date( idxEvents(i) );
        ShadowEvents.mag(iRow)              = Catalog.mag( idxEvents(i) );
        %
        ShadowEvents.magLatest(iRow)        = Catalog.mag( idxLatest );
        ShadowEvents.timeDiffLatest(iRow)   = Catalog.t(idxEvents(i)) - Catalog.t(idxLatest);
        ShadowEvents.distanceLatest(iRow)   = getDistance( Catalog.x(idxLatest), Catalog.y(idxLatest), ...
                                                           Catalog.x(idxEvents(i)), Catalog.y(idxEvents(i)), spaceUnit ) ...
                                              / estimate_rupSize(Catalog.mag(idxLatest), 'subduction', 'reverse', spaceUnit);
        %
        if ~(idxLatest==idxStrongest)
            ShadowEvents.magStrongest(iRow)     = Catalog.mag( idxStrongest );
            ShadowEvents.timeDiffStrongest(iRow)= Catalog.t(idxEvents(i)) - Catalog.t(idxStrongest);
            ShadowEvents.distanceStrongest(iRow)= getDistance( Catalog.x(idxStrongest), Catalog.y(idxStrongest), ...
                                                               Catalog.x(idxEvents(i)), Catalog.y(idxEvents(i)), spaceUnit ) ...
                                                  / estimate_rupSize(Catalog.mag(idxStrongest), 'subduction', 'reverse', spaceUnit);
        else
            ShadowEvents.idStrongest(iRow)      = NaN;
            ShadowEvents.magStrongest(iRow)     = NaN;
            ShadowEvents.timeDiffStrongest(iRow)= NaN;
            ShadowEvents.distanceStrongest(iRow)= NaN;
        end
        
    end
    
end
warning on