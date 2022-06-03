function doubletPercentages = group_doubletPerc4Sequences( tAnalysis, magnIntervals )

    %% Initialize auxiliary matrices and output matrix
    % Dimension of output matrix
    magnitudes          = tAnalysis.mainMw;
    nIntervals          = length(magnIntervals)-1;
    doubletPercentages  = zeros(1, nIntervals);
    rawPercentages      = sum(tAnalysis.percDoublets04_1(:,2:3), 2);
    
    % Loop over magnitude intervals
    for iInt = 1:nIntervals
        isInInterval = magnitudes >= magnIntervals(iInt) & magnitudes < magnIntervals(iInt+1);
        doubletPercentages(iInt) = mean( rawPercentages(isInInterval) );
    end
                                                         
end

