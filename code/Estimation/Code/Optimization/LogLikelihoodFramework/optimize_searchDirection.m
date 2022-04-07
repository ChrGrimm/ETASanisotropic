function [ searchDirection, ...
           hessMatrix, ...
           stopIteration ] = optimize_searchDirection( iInnerIter, ...
                                                         iOuterIter, ...
                                                         loglik_grad, ...
                                                         loglik_grad0, ...
                                                         hessMatrix, ...
                                                         dParam )
    
    %% Thresholds
    const1  = 1.0e-17;
    tau1    = 1.0e-6; 
    tau2    = 1.0e-6; 
    
    %% Initialize outputs
    searchDirection = zeros(length(dParam), 1);
    stopIteration   = false;
    
    if iInnerIter > 0 || iOuterIter > 1
           
          % Difference new to old gradient
          dLoglik_grad  = loglik_grad - loglik_grad0;
          wrk           = hessMatrix * dLoglik_grad; % H * dLoglik_grad

          sum1 = sum( wrk .* dLoglik_grad ); % dLoglik_grad' * H * dLoglik_grad'
          sum2 = sum( dParam .* dLoglik_grad );

          % Return criterion 1
          if sum1 <= const1 || sum2 <= const1
              disp('Return by sum1 <= const1 || sum2 <= const1')
              stopIteration = true;
              return;
          end

         if sum1 <= sum2
            % fletcher type correction
            stem = sum1/sum2 + 1; 
            for i = 1:length(dParam)
                for j = i:length(dParam)
                    hessMatrix(i,j) = hessMatrix(i,j) - (dParam(i) * wrk(j) + wrk(i) * dParam(j) - dParam(i) * dParam(j) * stem) / sum2;
                    hessMatrix(j,i) = hessMatrix(i,j);
                end
            end
         else
            for i = 1:length(dParam)
                for j = i:length(dParam)
                    hessMatrix(i,j) = hessMatrix(i,j) + dParam(i) * dParam(j) / sum2 - wrk(i) * wrk(j) / sum1;
                    hessMatrix(j,i) = hessMatrix(i,j);
                end
            end
         end

     end

%            ss = sum((hessMatrix*loglik_grad).^2); % H*loglik_grad
       searchDirection  = -hessMatrix * loglik_grad;
       sum1             = sum( searchDirection .* loglik_grad ); % -loglik_grad'*H*loglik_grad
       sum2             = sum( loglik_grad .* loglik_grad );
       sum2_sqrt        = sqrt( sum2);
       gtem             = abs( sum1 ) / sum2_sqrt;

       % Return criterion 2          
       if gtem <= tau1 && sum2_sqrt <= tau2              
           disp('Return by gtem <= tau1 && sum2_sqrt <= tau2')
           stopIteration = true;
           return;
       end

       if sum1 >= 0
           hessMatrix = diag(8*ones(8,1));
           searchDirection = -searchDirection;
       end
       
end
    
    