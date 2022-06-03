function [ magnIntervals_grouped, ...
           idxMagnIntervals2merge, ...
           strXtickLabels, ...
           idHistCat ] = define_plotSettings4doublets( strRegion, strVersion )

    % Define magnitude interval, which is used to group doublet percentages
    if strcmp(strRegion,'JPN') && strcmp(strVersion,'single plot')
        magnIntervals_grouped   = [5.9, 6.1, 6.3, 6.7, 8.7, 8.8, 9.7];
        idxMagnIntervals2merge  = [4,6];
        strXtickLabels          = {'[5.9,6.0]', '[6.1,6.2]','[6.3,6.6]','>=6.7','Tohoku'};
        idHistCat               = [1,3,5];
    elseif strcmp(strRegion,'JPN') && strcmp(strVersion,'shared plot')
        magnIntervals_grouped   = [5.9:0.2:7.1, 7.4:0.3:8.0, 9.7];
        idxMagnIntervals2merge  = [];
        strXtickLabels          = NaN;
        idHistCat               = [];
    elseif strcmp(strRegion,'CAL') && strcmp(strVersion,'single plot')
        magnIntervals_grouped   = [5.9, 6.1, 6.3, 6.7, 9.7];
        idxMagnIntervals2merge  = [];
        strXtickLabels          = {'[5.9,6.0]', '[6.1,6.2]','[6.3,6.6]','[6.7,M_{max}]'};
        idHistCat               = [2,4,6];
    elseif strcmp(strRegion,'CAL') && strcmp(strVersion,'shared plot')
        magnIntervals_grouped   = [5.9, 6.1, 6.3, 6.7, 9.7];
        idxMagnIntervals2merge  = [];
        strXtickLabels          = NaN;
        idHistCat               = [];
    end
    
%     if strcmp(strRegion,'JPN')
%         idHistCat = [1,3,5];
%     elseif strcmp(strRegion,'CAL')
%         idHistCat = [2,4,6];
%     end

end