function [ Analysis_detail, tAnalysis ]= analyze_sequences( Sequences, ...
                                                            DoubletCriteria, ...
                                                            fieldName, ...
                                                            spaceUnit, ...
                                                            magnThresh )

    %% Extract data from DoubletCriteria
    magnWindow          = round(DoubletCriteria.magnWindow,1);
    factorSpatExtent    = DoubletCriteria.factorSpatExtent;
                                                        
    %% Initialize output variables
    setOfSequences  = fieldName;
    clusterSize_1   = zeros(size(Sequences));
    clusterSize_2   = zeros(size(Sequences));
    meanTimeDiff    = zeros(size(Sequences));
    maxTimeDiff     = zeros(size(Sequences));
    meanDistance    = zeros(size(Sequences));
    maxDistance     = zeros(size(Sequences));
    bathLawMagn_1   = zeros(size(Sequences));
    bathLawMagn_2   = zeros(size(Sequences));
    dDoublet04_1    = zeros(length(Sequences),1);
    dDoublet04_2    = zeros(length(Sequences),1);
    
    %% Analyze sequences
    % Number of sequences
    nSequences      = length(Sequences);
    % Magnitude of sequence initiating mainshock
    mainshockMagn   = Sequences{1}.mag(1);
    
    % Loop over sequences
    for iSeq = 1:nSequences
        % Extract current sequence
        iSequence           = Sequences{iSeq};
        % Evaluate cluster size (not counting mainshock)
        clusterSize_2(iSeq) = size(iSequence,1)-1;

        %% Evaluate properties of sequences
        if clusterSize_2(iSeq)==0
            
            %% If no aftershocks were triggered, no properties are computed
            meanTimeDiff(iSeq)  = NaN;
            maxTimeDiff(iSeq)   = NaN;
            meanDistance(iSeq)  = NaN;
            maxDistance(iSeq)   = NaN;
            % Bath law
%             if mainshockMagn - magnThresh >= 1.5
%                 bathLawMagn_1(iSeq) = mainshockMagn - magnThresh + 0.1;
%                 bathLawMagn_2(iSeq) = mainshockMagn - magnThresh + 0.1;
%             else
                % set Bath law magnitude difference to NaN
                bathLawMagn_1(iSeq) = NaN;
                bathLawMagn_2(iSeq) = NaN;
%             end

        else       
            
            %% If aftershocks were triggered
            % Statistics of time differences
            tAftershocks        = iSequence.t(2:end);
            meanTimeDiff(iSeq)  = mean(tAftershocks);
            maxTimeDiff(iSeq)   = max(tAftershocks);
            
            % Statistics of distances
            iX                  = iSequence.x;
            iY                  = iSequence.y;
            distances           = getDistance( iX(2:end), iY(2:end), iX(1), iY(1), spaceUnit );
            meanDistance(iSeq)  = mean(distances);
            maxDistance(iSeq)   = max(distances);
            
            % Find events occurred within time-space window of mainshock
            isInSpace           = distances <= factorSpatExtent * iSequence.rupL(1);
            isInTime            = tAftershocks <= 365;
            isInSpaceAndTime    = isInSpace&isInTime;
            
            % Compute cluster size within time-space-window
            clusterSize_1(iSeq) = sum(isInSpaceAndTime);
            
            % Extract magnitudes of aftershocks
            aftershockMagn      = iSequence.mag(2:end);
            
            % Evaluate Bath law 
            if any(isInSpaceAndTime)
                bathLawMagn_1(iSeq)   = mainshockMagn - max(aftershockMagn(isInSpaceAndTime));
                bathLawMagn_2(iSeq)   = mainshockMagn - max(aftershockMagn);
            else
%                 if mainshockMagn - magnThresh >= 1.5
%                     bathLawMagn_1(iSeq) = mainshockMagn - magnThresh + 0.1;
%                     bathLawMagn_2(iSeq) = mainshockMagn - magnThresh + 0.1;
%                 else
                    bathLawMagn_1(iSeq) = NaN;
                    bathLawMagn_2(iSeq) = NaN;
%                 end
            end
            
            % Evaluate number of doublet partner events
            absMagnDifference   = round(abs(aftershockMagn-mainshockMagn),1);
            isInMagn2           = absMagnDifference <= magnWindow;
            
            dDoublet04_1(iSeq)  = sum( isInMagn2 & isInSpaceAndTime );
            dDoublet04_2(iSeq)  = sum( isInMagn2 );
            
        end
     
    end
       
    %% Store results in struct
    Analysis_detail = struct('setOfSequences', setOfSequences, 'clusterSize_1', clusterSize_1, 'clusterSize_2', clusterSize_2, ...
                             'meanTimeDiff', meanTimeDiff, 'maxTimeDiff', maxTimeDiff, ...
                             'meanDistance', meanDistance, 'maxDistance', maxDistance, ...
                             'bathLawMagn_1', bathLawMagn_1, 'bathLawMagn_2', bathLawMagn_2, ...
                             'dDoublet04_1', dDoublet04_1, 'dDoublet04_2', dDoublet04_2 );
    
    %% Compute average statistics
    tAnalysis                = table;
    tAnalysis.mainMw         = str2double(setOfSequences(end-1:end))/10;
    tAnalysis.clusterSize_1  = [quantile(clusterSize_1,0.05), mean(clusterSize_1), quantile(clusterSize_1,0.95)];
    tAnalysis.clusterSize_2  = [quantile(clusterSize_2,0.05), mean(clusterSize_2), quantile(clusterSize_2,0.95)];
    tAnalysis.meanTimeDiff   = [quantile(meanTimeDiff,0.05), mean(meanTimeDiff,'omitnan'), quantile(meanTimeDiff, 0.95)];
    tAnalysis.maxTimeDiff    = [quantile(maxTimeDiff,0.05), mean(maxTimeDiff,'omitnan'), quantile(maxTimeDiff, 0.95)];
    tAnalysis.meanDistance   = [quantile(meanDistance,0.05), mean(meanDistance,'omitnan'), quantile(meanDistance, 0.95)];
    tAnalysis.maxDistance    = [quantile(maxDistance,0.05), mean(maxDistance,'omitnan'), quantile(maxDistance, 0.95)];
    bathLawMagn_1            = bathLawMagn_1(~isnan(bathLawMagn_1));
    bathLawMagn_2            = bathLawMagn_2(~isnan(bathLawMagn_2));
    tAnalysis.bathLawMagn_1  = [quantile(bathLawMagn_1,0.05), mean(bathLawMagn_1), quantile(bathLawMagn_1, 0.95)];
    tAnalysis.bathLawMagn_2  = [quantile(bathLawMagn_2,0.05), mean(bathLawMagn_2), quantile(bathLawMagn_2, 0.95)];
    
    tAnalysis.percDoublets04_1   = [sum(dDoublet04_1==0), sum(dDoublet04_1==1), sum(dDoublet04_1>1)] / nSequences;
    tAnalysis.percDoublets04_2   = [sum(dDoublet04_2==0), sum(dDoublet04_2==1), sum(dDoublet04_2>1)] / nSequences;

end