function plot_spatialDistr( modelDesign, paramETAS, M_c, strFolder )
% 
% This function plots the spatial distribution for two exemplary magnitudes
% resulting from an ETAS parameter estimation beforehand.
%
% Call: plot_spatialDistr( paramETAS, modelDesign, M_c )
%
% INPUT:
%
% modelDesign:  string variable         controlling the chosen model design
% paramETAS:    numeric vector (8x1)    of estimated ETAS parameters
% M_c:           scalar                  (lower) magnitude threshold
% strFolder:    string variable         name of folder describing scenario 

    % Extract relevant ETAS parameters
    D       = paramETAS(6);
    gamma   = paramETAS(7);
    q       = paramETAS(8);
    
    % Sample of two exemplary magnitudes (reduced by threshold M_c)
    magn = [6.5, 8.5] - M_c;
    
    % Precompute f_inn
    % !!! Outdated inputs !!!
    [ ~, f_inn, ~ ] = set_modelFunctions( modelDesign, sqrt(paramETAS), magn );
    
    % Prepare grid
    xyGrid          = linspace(-1, 1, 100);
    [meshX, meshY]  = meshgrid(xyGrid, xyGrid);
    
    if strcmp(modelDesign, 'standard-iso')
        
        density1 = (q-1)./(pi*D*exp(gamma*magn(1))) * f_inn(meshX.^2+meshY.^2, magn(1)).^(-q);
        density2 = (q-1)./(pi*D*exp(gamma*magn(2))) * f_inn(meshX.^2+meshY.^2, magn(2)).^(-q);
        plot3_nested(meshX(:), meshY(:), density1, magn(1), strFolder)
        plot3_nested(meshX(:), meshY(:), density2, magn(2), strFolder)
        
    elseif strcmp(modelDesign, 'standard-iso 2')
        
        density1 = (q-1)./(D*exp(gamma*magn(1))) * f_inn(meshX.^2+meshY.^2, magn(1)).^(-q);
        density2 = (q-1)./(D*exp(gamma*magn(2))) * f_inn(meshX.^2+meshY.^2, magn(2)).^(-q);
        plot3_nested(meshX(:), meshY(:), density1, magn(1), strFolder)
        plot3_nested(meshX(:), meshY(:), density2, magn(2), strFolder)
        
    elseif strcmp(modelDesign, 'aniso')
        
        factor  = (q-1)/(D*exp(gamma*magn));
        rupL    = get_rupLength( magn+4.5, 2 );
        r1      = get_distanceLine( [-rupL(1)/2, 0], [rupL(1), 0], [meshX, meshY], 0 );
        density1 = (q-1)/(D*exp(gamma*magn)) * f_inn(r1, r1.^2, rupL(1), magn(1)).^(-q);
        r2      = get_distanceLine( [-rupL(2)/2, 0], [rupL(2), 0], [meshX, meshY], 0 );
        density2 = (q-1)/(D*exp(gamma*magn)) * f_inn(r2, r2.^2, rupL(2), magn(2)).^(-q);
        plot3_nested(meshX(:), meshY(:), density1, magn(1), strFolder)
        plot3_nested(meshX(:), meshY(:), density2, magn(2), strFolder)
        
        
    end
    
    function plot3_nested(x,y,z,magn, strFolder)
        figure;
        plot3(x,y,z(:))
        xlabel('Space (in degrees)')
        zlabel('PDF of offspring occurrence')
        title([strFolder, ': Spatial PDF for Magn = ', num2str(M_c+magn)])
    end

end

