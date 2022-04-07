function check_weightsSpatDistr( weight, iEvent, strSegmType )
% This function checks if the segment weights of the precomputed isotropic/anisotropic spatial
% distributions of the ETAS model sum up to the required total weight.
%
% Call: check_weightsSpatDistr( weight, iEvent, strSegmType )
%
% INPUTS:
%
% - weight:         numeric vector (Nx1)    spatial weight per radial or box segment
% - iEvent:         scalar                  event ID (for tracking purpose in error message)
% - strSegmType:    string                  description of segment type (e.g. 'radial' or 'anisotropic')

    warning off backtrace
    eps = 10^(-6);
    
    if strcmp(strSegmType, 'radial (aniso mode)')
        
        % Radial weights must sum up to either 0 or 0.5
        if abs( sum(weight)-0.5 ) > sqrt(eps) && abs( sum(weight)-0 ) > sqrt(eps)
            
            error(['Event ', num2str(iEvent), ' has radial weights summing up to ', num2str(sum(weight)), '. Should be 0 or 0.5!'])
            
        elseif abs( sum(weight)-0.5 ) > eps && abs( sum(weight)-0 ) > eps
            
            warning(['ERROR? Event ', num2str(iEvent), ' has radial weights summing up to ', num2str(sum(weight)), '. Should be 0 or 0.5!'])

        end
        
    elseif strcmp(strSegmType, 'radial (iso mode)')
    
        % Radial weights must sum up to either 0 or 1
        if abs( sum(weight)-0.5 ) < eps
            
            warning(['Event ', num2str(iEvent), ' has radial weights summing up to 0.5, i.e. event location must lie on polygon edge!'])
            
        elseif abs( sum(weight)-1 ) > eps && abs( sum(weight)-0 ) > eps
            
            error(['Event ', num2str(iEvent), ' has radial weights not summing up to 0 or 1!'])

        end
        
    end

end