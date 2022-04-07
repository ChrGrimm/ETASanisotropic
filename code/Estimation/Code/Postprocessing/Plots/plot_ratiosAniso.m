i = 1;
iMag = magnitudes(i)-M_c;
sampleR = 0:0.01:10;

% Rupture length
rupL = get_rupLength(iMag+M_c, typeRupL) /111;

% Compute ratio sample
F_cdf       = @(r) 1 - f_inn(r, r.^2, rupL, iMag).^(1-q);
f_anisoDens = @(r) (q-1)/(D*exp(gamma*iMag)) * (2*rupL) .* f_inn(r, r.^2, rupL, iMag).^(-q);
f_isoDens   = @(r) (q-1)/(D*exp(gamma*iMag)) * (2*pi*r) .* f_inn(r, r.^2, rupL, iMag).^(-q);
f_ges       = @(r) f_isoDens(r)+f_anisoDens(r);
for i=1:length(sampleR)
    F_ges(i) = integral(f_ges, 0, sampleR(i));
end
max(abs(F_cdf(sampleR)-F_ges))

% Ratio sample for derivative by D (same as for gamma!)
F_dD        = @(r) (1-q)/D .* f_inn(r, r.^2, rupL, iMag).^(1-q) .* (1 - f_inn(r, r.^2, rupL, iMag).^(-1));
f_aniso_dD  = @(r) f_anisoDens(r)/D .* (q * (1-f_inn(r,r.^2,rupL,iMag).^(-1)) - 1);
f_iso_dD    = @(r) f_isoDens(r)/D .* (q * (1-f_inn(r,r.^2,rupL,iMag).^(-1)) - 1);
for i=1:length(sampleR)
    F_iso_dD(i)     = integral(f_iso_dD, 0, sampleR(i)); 
    F_aniso_dD(i)   = integral(f_aniso_dD, 0, sampleR(i));
end
max(abs(F_dD(sampleR)-(F_aniso_dD+F_iso_dD)))

figure;
plot(sampleR, F_dD(sampleR))
hold on
plot(sampleR, F_aniso_dD)
plot(sampleR, F_iso_dD)
legend('F dD', 'F aniso dD', 'F iso dD', 'Location', 'Northeast')

%%
% Rupture length
rupL = rupL/2;

% Compute ratio sample
F_cdf       = @(r) 1 - f_inn(r, r.^2, rupL, iMag).^(1-q);
f_anisoDens = @(r) (q-1)/(D*exp(gamma*iMag)) * (2*rupL) .* f_inn(r, r.^2, rupL, iMag).^(-q);
f_isoDens   = @(r) (q-1)/(D*exp(gamma*iMag)) * (2*pi*r) .* f_inn(r, r.^2, rupL, iMag).^(-q);
f_ges       = @(r) f_isoDens(r)+f_anisoDens(r);
for i=1:length(sampleR)
    F_ges(i) = integral(f_ges, 0, sampleR(i));
end
max(abs(F_cdf(sampleR)-F_ges))

% Ratio sample for derivative by D (same as for gamma!)
F_dD        = @(r) (1-q)/D .* f_inn(r, r.^2, rupL, iMag).^(1-q) .* (1 - f_inn(r, r.^2, rupL, iMag).^(-1));
f_aniso_dD  = @(r) f_anisoDens(r)/D .* (q * (1-f_inn(r,r.^2,rupL,iMag).^(-1)) - 1);
f_iso_dD    = @(r) f_isoDens(r)/D .* (q * (1-f_inn(r,r.^2,rupL,iMag).^(-1)) - 1);
for i=1:length(sampleR)
    F_iso_dD(i)     = integral(f_iso_dD, 0, sampleR(i)); 
    F_aniso_dD(i)   = integral(f_aniso_dD, 0, sampleR(i));
end
max(abs(F_dD(sampleR)-(F_aniso_dD+F_iso_dD)))

figure;
plot(sampleR, F_dD(sampleR))
hold on
plot(sampleR, F_aniso_dD)
plot(sampleR, F_iso_dD)
legend('F dD', 'F aniso dD', 'F iso dD', 'Location', 'Northeast')

%%
% Rupture length
rupL = rupL/100;

% Compute ratio sample
F_cdf       = @(r) 1 - f_inn(r, r.^2, rupL, iMag).^(1-q);
f_anisoDens = @(r) (q-1)/(D*exp(gamma*iMag)) * (2*rupL) .* f_inn(r, r.^2, rupL, iMag).^(-q);
f_isoDens   = @(r) (q-1)/(D*exp(gamma*iMag)) * (2*pi*r) .* f_inn(r, r.^2, rupL, iMag).^(-q);
f_ges       = @(r) f_isoDens(r)+f_anisoDens(r);
for i=1:length(sampleR)
    F_ges(i) = integral(f_ges, 0, sampleR(i));
end
max(abs(F_cdf(sampleR)-F_ges))

% Ratio sample for derivative by D (same as for gamma!)
F_dD        = @(r) (1-q)/D .* f_inn(r, r.^2, rupL, iMag).^(1-q) .* (1 - f_inn(r, r.^2, rupL, iMag).^(-1));
f_aniso_dD  = @(r) f_anisoDens(r)/D .* (q * (1-f_inn(r,r.^2,rupL,iMag).^(-1)) - 1);
f_iso_dD    = @(r) f_isoDens(r)/D .* (q * (1-f_inn(r,r.^2,rupL,iMag).^(-1)) - 1);
for i=1:length(sampleR)
    F_iso_dD(i)     = integral(f_iso_dD, 0, sampleR(i)); 
    F_aniso_dD(i)   = integral(f_aniso_dD, 0, sampleR(i));
end
max(abs(F_dD(sampleR)-(F_aniso_dD+F_iso_dD)))

figure;
plot(sampleR, F_dD(sampleR))
hold on
plot(sampleR, F_aniso_dD)
plot(sampleR, F_iso_dD)
legend('F dD', 'F aniso dD', 'F iso dD', 'Location', 'Northeast')