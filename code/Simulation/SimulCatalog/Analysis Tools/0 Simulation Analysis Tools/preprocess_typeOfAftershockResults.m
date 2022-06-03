function [ centerMagnitudes, ...
           resultsBath, ...
           resultsDoublet ] = preprocess_typeOfAftershockResults( tAnalysis_synthCat )

    %% Initializations
    % Auxiliary variables and matrices
    eventRelationBath       = zeros(5000000,2);
    eventRelationDoublets   = zeros(500000,2);
    iCounterBath            = 0;
    iCounterDoublets        = 0;
    nRealiz                 = size(tAnalysis_synthCat,1);
    % Output vectors and matrices
    uMagnitudes     = [5.5:0.1:7.5, 7.7:0.2:8.1, 8.7, 8.8]; % , max(eventRelationBath(:,1))+0.1];
    centerMagnitudes = (uMagnitudes(1:end-1)+uMagnitudes(2:end))/2;
    resultsBath     = zeros(length(uMagnitudes)-1,3);
    resultsDoublet  = zeros(length(uMagnitudes)-1,3);
    
    %% Compile bath law and doublet event relations in single matrices
    for iRealiz = 1:nRealiz
        nBath                                                   = size(tAnalysis_synthCat.bathLawType{iRealiz},1);
        nDoublets                                               = size(tAnalysis_synthCat.doubletType{iRealiz},1);
        eventRelationBath(iCounterBath+(1:nBath),:)             = tAnalysis_synthCat.bathLawType{iRealiz};
        eventRelationDoublets(iCounterDoublets+(1:nDoublets),:) = tAnalysis_synthCat.doubletType{iRealiz};
        iCounterBath                                            = iCounterBath + nBath;
        iCounterDoublets                                        = iCounterDoublets + nDoublets;
    end
    % Delete empty, pre-initialized rows
    eventRelationBath(iCounterBath+1:end,:)         = [];
    eventRelationDoublets(iCounterDoublets+1:end,:) = [];
   
    %% Compute percentages per unique magnitude step    
    for i=1:length(uMagnitudes)-1
        isThisMagn          = eventRelationBath(:,1)>=uMagnitudes(i) & eventRelationBath(:,1)<uMagnitudes(i+1);
%         resultsBath(i,:)    = [ sum(eventRelationBath(isThisMagn,2)==0), ...
%                                 sum(ismember(eventRelationBath(isThisMagn,2),[1,2])), ...
%                                 sum(eventRelationBath(isThisMagn,2)==3) ] / sum(isThisMagn);
        resultsBath(i,:)    = [ sum(eventRelationBath(isThisMagn,2)==0), ...
                                sum(eventRelationBath(isThisMagn,2)==1), ...
                                sum(eventRelationBath(isThisMagn,2)>=2) ] / sum(isThisMagn);
        isThisMagn          = eventRelationDoublets(:,1)>=uMagnitudes(i) & eventRelationDoublets(:,1)<uMagnitudes(i+1);
%         resultsDoublets(i,:)= [ sum(eventRelationDoublets(isThisMagn,2)==0), ...
%                                 sum(ismember(eventRelationDoublets(isThisMagn,2),[1,2])), ...
%                                 sum(eventRelationDoublets(isThisMagn,2)==3) ] / sum(isThisMagn);
        resultsDoublet(i,:) = [ sum(eventRelationDoublets(isThisMagn,2)==0), ...
                                sum(eventRelationDoublets(isThisMagn,2)==1), ...
                                sum(eventRelationDoublets(isThisMagn,2)>=2) ] / sum(isThisMagn);
    end
    
end