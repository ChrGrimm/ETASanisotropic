function [ multipletDegree, ...
           idxEvents, ...
           multipletType, ...
           EventRelations_bath, ...
           magnitudePairs ] = remove_nanRows4doubletAnalysis( multipletDegree, ...
                                                              idxEvents, ...
                                                              multipletType, ...
                                                              EventRelations_bath, ...
                                                              magnitudePairs )
                                                          
    % Remove NaN rows from multipletDegree and idxEvents
    isnanDegrees                        = isnan(multipletDegree(:,1));
    multipletDegree(isnanDegrees,:)     = [];
    idxEvents(isnanDegrees)             = [];
    % Remove NaN rows from multipletType
    isNanData                           = isnan(multipletType(:,1));
    multipletType(isNanData,:)          = [];
    % Remove NaN rows from EventRelations_bath and magnitudePairs
    isNanData                           = isnan(EventRelations_bath(:,1));
    EventRelations_bath(isNanData,:)    = [];
    magnitudePairs(isNanData,:)         = [];                                                      
                                                          
end