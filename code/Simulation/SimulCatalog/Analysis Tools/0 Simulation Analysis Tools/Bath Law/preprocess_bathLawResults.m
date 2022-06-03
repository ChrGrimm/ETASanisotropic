function [ bathLawSeq_1, ...
           bathLawSeq_2, ...
           bathLawCat ] = preprocess_bathLawResults( tAnalysis_synthSeq, ...
                                                     tAnalysis_synthCat )
     
    %% Filter data by magnitude
    isMagn          = tAnalysis_synthSeq.mainMw >= 5.5;
    magnitudes      = tAnalysis_synthSeq.mainMw(isMagn);
    bathLawSeq_1    = tAnalysis_synthSeq.bathLawMagn_1(isMagn, :);
    bathLawSeq_2    = tAnalysis_synthSeq.bathLawMagn_2(isMagn, :);
    
    %% Compute bath law data for synthetic catalogs
    bathLawData     = tAnalysis_synthCat.bathLawMagnDiff;
    bathLawAll      = zeros(2*10000*size(bathLawData{1},1), 3);
    c = 0;
    for j = 1:size(bathLawData,1)
        jBathLawData                    = bathLawData{j};
        lengthData                      = size(jBathLawData,1);
        bathLawAll(c+1:c+lengthData, :) = jBathLawData;
        c                               = c + lengthData;
    end
    bathLawAll( bathLawAll(:,2)<0, : )  = [];
    bathLawCat = zeros(length(magnitudes),3);
    for j = 1:length(magnitudes)
         jDiffData      = bathLawAll( bathLawAll(:,1) == magnitudes(j), :);
         jDiffData      = repelem(jDiffData(:,2), round(jDiffData(:,3)));
         bathLawCat(j,:) = [ quantile(jDiffData,0.1), mean(jDiffData), quantile(jDiffData,0.9) ];
    end                                          
                                                 
end