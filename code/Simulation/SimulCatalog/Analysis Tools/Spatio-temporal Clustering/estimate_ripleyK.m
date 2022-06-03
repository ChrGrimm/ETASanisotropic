function K = estimate_ripleyK( x, y, A )
  
    idx_rand= randperm(length(x));
    N       = min(2500, length(x));
    x       = x(idx_rand(1:N));
    y       = y(idx_rand(1:N));
    gamma   = A*N^-2;
    dists   = getDistance( x, y, x', y', 'degree' );
    dists   = dists(:);
    h       = 10.^(-1:.5:3) / 111; % in degrees
    K       = gamma * ( sum(dists<=h, 1)-N )'; 

end