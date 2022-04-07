function [ strikes, ...
            epiPos, ...
            nAftershocks ] = find_maxTrigDensity( forwardRate_sum, ...
                                                  nAftershocks_matrix, ...
                                                  sampleStrike, ...
                                                  sampleEpiPos, ...
                                                  isDuplicate )

    %% Initializations
    nEvents     = size(forwardRate_sum,1);
    strikes     = zeros(nEvents, 1);
    epiPos      = zeros(nEvents, 1);
    nAftershocks= zeros(nEvents, 1);
    forwardRate = zeros(nEvents, 1);
    
    %% Find maximum
    for iEvent = 1:nEvents
        iForwardRates = squeeze( forwardRate_sum(iEvent,:,:) );
        
        % Maximize forward rates per strike angle over epicenter locations
        maxForwardRates = max(iForwardRates, [],2);
        
        % Identify local maxima in forward rates as a function of strike angle
        idxLocalMaxima = find( islocalmax(maxForwardRates) ); % imregionalmax
        if maxForwardRates(1) > max( maxForwardRates([2,end]) )
            idxLocalMaxima = [idxLocalMaxima; 1];
        end
        if maxForwardRates(end) > max( maxForwardRates([1,end-1]) )
            idxLocalMaxima = [idxLocalMaxima; size(maxForwardRates,1)];
        end
        if isempty(idxLocalMaxima); idxLocalMaxima = 1; end
        [ ~, idxSort ] = sort( maxForwardRates(idxLocalMaxima), 'descend' );
        
        % Store strike angle and epicenter location
        if isDuplicate(iEvent)
            idxStrikeAngle = idxLocalMaxima(idxSort(2));
        else
            idxStrikeAngle = idxLocalMaxima(idxSort(1));
        end
        strikes(iEvent)     = sampleStrike( iEvent, idxStrikeAngle );
        forwardRate(iEvent) = maxForwardRates( idxStrikeAngle );
        [~, idxEpiPosition] = max( iForwardRates(idxStrikeAngle,:) );
        if idxEpiPosition==1
            idxSameRate     = find(iForwardRates(idxStrikeAngle,:)==iForwardRates(idxStrikeAngle,1));
            diffToCentered  = abs( sampleEpiPos(idxSameRate)-0.5 );
            idxEpiPosition  = idxSameRate( diffToCentered==min(diffToCentered) );
            idxEpiPosition  = idxEpiPosition(1);
        end
        
        epiPos(iEvent)      = sampleEpiPos( idxEpiPosition );
        nAftershocks(iEvent)= nAftershocks_matrix(iEvent, idxStrikeAngle, idxEpiPosition);
        
        %% Plotting%      
        makePlots = false;
%         makePlots = true;

        if makePlots
            % Scatter Plot
            if iEvent==2
                figure
                plot( sampleStrike(iEvent,:), iForwardRates, 'Color', 'r', 'LineWidth', 0.1 )
                hold on
                text( strikes(1)+10, 230, ['strike   = ', num2str(strikes(1)), '°'] )
                text( strikes(1)+10, 224, ['epiPos = ', num2str(epiPos(1))] )
                text( strikes(2)+10, 230, ['strike   = ', num2str(strikes(2)), '°'] )
                text( strikes(2)+10, 224, ['epiPos = ', num2str(epiPos(2))] )
                make_titleLeftCorner('(a)')
                
            elseif iEvent == 3
                figure
                plot( sampleStrike(iEvent,:), iForwardRates, 'Color', 'b', 'LineWidth', 0.1 )
                hold on
                text( strikes(3)-45, 925, ['strike   = ', num2str(strikes(3)), '°'] )
                text( strikes(3)-45, 890, ['epiPos = ', num2str(epiPos(3))] )
                make_titleLeftCorner('(b)')
            end
            
            xlabel('Strike angle (°)')
            ylabel('Forward trigger contribution')
%             legend('curves by epicenter location')

        end
    end

end

%% Old code pieces
% % Find index of maximum rate in vectorized matrix
% [~, idxMax]         = max( forwardRate_sum(:) );
% % Re-convert index of maximum to matrix indeces
% [ iStrike, iEpiPos ]= ind2sub( size(forwardRate_sum), idxMax );
% strikes(iEvent)     = sampleStrike( iEvent, iStrike );
% epiPos(iEvent)      = sampleEpiPos( iEpiPos );

% id2plot = [22, 671]; %[219543, 220235];
% isEvent = ismember( Catalog.id(idxAniso), id2plot );
% plot_rateVsStrike( sampleStrike(isEvent,:), ...
%                            sampleEpiPos, ...
%                            forwRate_sum(isEvent,:,:), ...
%                            forwRate_logl(isEvent,:,:), ...
%                            '-' )
