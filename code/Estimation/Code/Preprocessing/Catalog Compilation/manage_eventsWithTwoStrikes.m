function Catalog = manage_eventsWithTwoStrikes( Catalog, evIDs_twoStrikes, version )

    if any( ismember(Catalog.id, evIDs_twoStrikes) )
        if ismember(version, {'estimation', 'simulation'})
            idx                     = find( ismember(Catalog.id, evIDs_twoStrikes) );
                                        
            if ~ all( strcmp(Catalog.typeKernel(idx), 'aniso') )
                error('An event that should be modeled with two strikes, is not anisotropic. Check eventID and SpaceSettings.anisoFromMw.')
            end
                
            for i = length(idx):-1:1
                line2duplicate              = Catalog( idx(i), : );
                Catalog                     = [ Catalog(1:idx(i),:); ...
                                                line2duplicate; ...
                                                Catalog(idx(i)+1:size(Catalog,1),:) ];
                Catalog.isDupl(idx(i)+1)    = true;
                Catalog.evWeight(idx(i))    = 0.5;
                Catalog.evWeight(idx(i)+1)  = 0.5;
            end


        elseif strcmp(version, 'merge')
            idx                         = find( Catalog.id == evIDs_twoStrikes ); 
            Catalog.backgrRate(idx(1))  = sum( Catalog.backgrRate(idx) );
            Catalog.lambda(idx(1))      = sum( Catalog.lambda(idx) );
            Catalog.evWeight(idx(1))    = 1;
%             Catalog.spatIntegr(idx(1))  = sum( Catalog.spatIntegr(idx) .* Catalog.evWeight(idx) );
            Catalog(idx(2),:)           = [];

        end

    end
end

%% Old code pieces
%     elseif strcmp(version, 'simulation')
%         idx                     = find( Catalog.id == evIDs_twoStrikes );
%         for i = length(idx):-1:1
%             line2duplicate          = Catalog( idx(i), : );
%             Catalog                 = [ Catalog(1:idx(i),:); ...
%                                         line2duplicate; ...
%                                         Catalog(idx(i)+1:size(Catalog,1),:) ];
%             Catalog.id(idx(i)+1)       = evIDs_twoStrikes + 0.1;
%             Catalog.strike(idx)     = SpaceSettings.duplicate_strikes(1);
%             Catalog.strike(idx+1)   = SpaceSettings.duplicate_strikes(2);
%             Catalog.epiPos(idx)     = SpaceSettings.duplicate_epiPos(1);
%             Catalog.epiPos(idx+1)   = SpaceSettings.duplicate_epiPos(2);
%             Catalog.evWeight(idx)   = 0.5;
%             Catalog.evWeight(idx+1) = 0.5;
%         end
