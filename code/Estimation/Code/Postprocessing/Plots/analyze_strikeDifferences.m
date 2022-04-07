function analyze_strikeDifferences( tempCat, strike )

    isMagn5 = tempCat.mag>=5.0;
    tempCat = tempCat(isMagn5,:);
    strike  = strike(isMagn5);

    StrikeCompare = [tempCat.str1, tempCat.str2, strike];
    StrikeCompare(StrikeCompare>=180) = StrikeCompare(StrikeCompare>=180)-180; 
    StrikeCompare(:,4) = min([abs(StrikeCompare(:,3)-StrikeCompare(:,1)), abs(StrikeCompare(:,3)-(StrikeCompare(:,1)+180)), abs(StrikeCompare(:,3)-(StrikeCompare(:,1)-180))], [], 2);
    StrikeCompare(:,5) = min([abs(StrikeCompare(:,3)-StrikeCompare(:,2)), abs(StrikeCompare(:,3)-(StrikeCompare(:,2)+180)), abs(StrikeCompare(:,3)-(StrikeCompare(:,2)-180))], [], 2);

    if min(tempCat.mag)==4.0
        magnRanges = [5.0, 5.5, 7.0, 10];
    else
        magnRanges = [5.0, 6.0, 7.5];
    end
    
    figure
    eps = 10^-60;
    for iRange=1:length(magnRanges)-1
        
        idxEvents = tempCat.mag >= magnRanges(iRange) & tempCat.mag < magnRanges(iRange+1);
        subplot(2,2,iRange)
        histogram( min(StrikeCompare(idxEvents,4), StrikeCompare(idxEvents,5)), -0.5:1:89.5, 'Normalization', 'probability' )
        hold on
        ksdensity( 0.1+min(StrikeCompare(idxEvents,4), StrikeCompare(idxEvents,5)), 0:90, 'support', 'positive' )
        tmpDensity = ksdensity( 0.1+min(StrikeCompare(idxEvents,4), StrikeCompare(idxEvents,5)), 0:90, 'support', 'positive' );
        plot(0:90, tmpDensity, 'Color', 'r', 'LineWidth', 2)
        title(['Kernel density: ', num2str(magnRanges(iRange)), ' <= Mw < ', num2str(magnRanges(iRange+1))])
        xlim([0, 90])
        
    end
    
    figure
    histogram( min(StrikeCompare(:,4), StrikeCompare(:,5)), -0.5:1:89.5, 'Normalization', 'probability' )
    hold on
    tmpDensity = ksdensity( 0.1+min(StrikeCompare(idxEvents,4), StrikeCompare(idxEvents,5)), 0:90, 'support', 'positive' );
    plot(0:90, tmpDensity, 'Color', 'r', 'LineWidth', 2)
    title('Kernel density: All magnitudes')
    xlim([0, 90])
              
%     [M,D]=meshgrid(unique(Catalog.mag), unique(min(StrikeCompare(:,4), StrikeCompare(:,5))));
%     for i=1:length(M(:))
%         nEvents(i) = sum(Catalog.mag==M(i) & min(StrikeCompare(:,4), StrikeCompare(:,5))==D(i));
%     end
% 
%     isOne = nEvents==1;
%     isTwo = nEvents==2;
%     isThreeToFive = nEvents>2 & nEvents<6;
%     isSixToTen = nEvents>5 & nEvents<11;
%     isAbove10 = nEvents>10;
% 
%     figure
%     scatter(M(isOne),D(isOne),5,'b','filled');
%     hold on
%     scatter(M(isTwo),D(isTwo),5,'r','filled');
%     scatter(M(isThreeToFive),D(isThreeToFive),5,'y','filled');
%     scatter(M(isSixToTen),D(isSixToTen),5,'g','filled');
%     scatter(M(isAbove10),D(isAbove10),5,'m','filled');
%     legend('1','2','3-5','6-10','>10')
    
end