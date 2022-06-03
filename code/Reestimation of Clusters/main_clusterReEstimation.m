function main_clusterReEstimation( minSizeCluster )

    % Load results from ETAS fit over the entire catalog
    Inputs = load('Inputs_estimation.mat');
    load('Results_estimation.mat', 'Catalog', 'ClustersETAS', 'paramETAS')
    
    idxLargeClusters = find( ClustersETAS.size >= minSizeCluster );
    
    %% Loop over large clusters
    for iCluster = idxLargeClusters'
        strFolderName = ['Cluster_', num2str( ClustersETAS.clusterID(iCluster) )];
        mkdir(strFolderName)
        addpath(strFolderName)
        
        % Create cluster specific inputs
        cd(strFolderName)
        modify_inputs4clusterReEstimation( Catalog, ClustersETAS(iCluster,:) , Inputs, paramETAS );
        
        disp(['Re-estimate Cluster ', num2str(iCluster)])
        
        main_etasEstimation( true )
        
        % Remove path
%         rmpath(strFolderName)
        cd('..')
        
    end

end