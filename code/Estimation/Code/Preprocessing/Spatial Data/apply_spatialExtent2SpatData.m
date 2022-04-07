function [r, r_squ, r0, r0_squ, segmFactor, isStartP, isEndP] ...
       = apply_spatialExtent2SpatData(spatialExtent, r, r0, segmFactor, isStartP, isEndP)
   
    % Boolean, which distances r0 are within spatial extent
    isInExtent      = r0 <= spatialExtent;

    % Update r0, r0_squ and segmFactor:
    % - Aggregate all segments outside of spatial extent to first value
    % - Reduce rest of the list to segments located within spatial extent only
    r0              = [ spatialExtent; r0(isInExtent) ];
    r0_squ          = r0.^2;
    segmFactor      = [ sum(segmFactor(~isInExtent)); segmFactor(isInExtent) ];
    % Update r, r_squ, isStartP and isEndP:
    % - Aggregate all segments outside of spatial extent to first value
    % - Reduce rest of the list to segments located within spatial extent only
    idxStartP       = find(isStartP);
    idxEndP         = find(isEndP);
    idxStartP       = idxStartP(isInExtent);
    idxEndP         = idxEndP(isInExtent);
    idxMerge        = unique([idxStartP;idxEndP]);
    % Take min here, since possibly r0 <= spatialExtent, but belonging r >
    % spatialExtent
    r               = [spatialExtent; min(spatialExtent, r(idxMerge))]; 
    r_squ           = r.^2;
    isStartP        = [true; ismember(idxMerge, idxStartP)];
    isEndP          = [true; ismember(idxMerge, idxEndP)];
    
end