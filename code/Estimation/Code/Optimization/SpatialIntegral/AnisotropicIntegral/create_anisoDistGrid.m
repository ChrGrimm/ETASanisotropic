function [ distGrid_integral, ...
           distGrid_gradiants ] = create_anisoDistGrid( peak_integr, ...
                                                        peak_grad, ...
                                                        nDistGrid_integr, ...
                                                        nDistGrid_grad, ...
                                                        R, ...
                                                        inclGradiants )

    %% Compute grid of spatial integral
%     distGrid_integral = [10^(-30), 10^(-8), 0.01:0.01:0.99, 1:0.1:2.9, 3:1:20, 1000];
    if peak_integr >= R
        % Case: Peak is beyond spatial extent
        distGrid_integral = unique([ 10^(-30), linspace(10^-8, R, nDistGrid_integr-1)]);
    else
        % Case: Peak is before spatial extent
        distGrid_integral = unique([ 10^(-30), linspace(10^-8, peak_integr, round(2/3*nDistGrid_integr)), ...
                                        linspace(peak_integr, R, nDistGrid_integr-round(2/3*nDistGrid_integr) )]);
    end
    
    %% Compute grid for spatial gradiants
%     gridR_dF(:,i) = [10^(-30), 10^(-8), 0.01:0.01:0.99, 1:0.1:2.9, 3:1:20, 1000];
    if inclGradiants
        if 2*peak_grad >= R
            % Case: Peak is beyond half of spatial extent
            distGrid_gradiants = unique([ 10^(-30), linspace(10^-8, R, nDistGrid_grad-1)]);
        else
            % Case: Peak is before half of spatial extent
            distGrid_gradiants = unique([ 10^(-30), linspace(10^-8, 2*peak_grad, round(2/3*nDistGrid_grad)), ...
                                           linspace(2*peak_grad, (2*peak_grad+R)/2, round(1/5*nDistGrid_grad)), ...
                                           linspace((2*peak_grad+R)/2, R, nDistGrid_grad-round(2/3*nDistGrid_grad)-round(1/5*nDistGrid_grad)+1)]);
        end
    else
        distGrid_gradiants = -1;
    end
    
end
