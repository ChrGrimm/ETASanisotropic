function magn = sample_magnitudesByGR( nEvents, Mc, Mmax, beta )

    % Sample magnitudes from Gutenberg-Richter
    magn        = Mc + round(exprnd(1/beta, nEvents, 1), 1);

    % Re-sample magnitudes larger than upper threshold
    while any(magn > Mmax)
        isLarger        = magn > Mmax;
        magn(isLarger)  = Mc + round(exprnd(1/beta, sum(isLarger), 1), 1); %- log(1-randN) / beta;
    end

end
