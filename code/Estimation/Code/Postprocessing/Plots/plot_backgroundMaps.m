figure;
scatter(Catalog.xx(notFlag2), Catalog.yy(notFlag2), 5, log(bwd), 'filled'); colorbar
idxNot2 = find(notFlag2); idx = idxNot2(bwd==0.05); hold on; scatter(Catalog.xx(idx), Catalog.yy(idx), 5, 'red', 'filled'); colorbar
title('Log(bandwidths) per event')
xlabel('Centered longitude (Japan)')
ylabel('Centered latitude (Japan)')
legend('Bandwidth>0.05', 'Bandwidth=0.05', 'Location', 'Northwest')
axis([-8 8 -10 10])

figure;
scatter(Catalog.xx(notFlag2), Catalog.yy(notFlag2), 5, log(sum(dGauss,2)), 'filled'); colorbar
title('Log(Sum of distributed background probability density) per event')
xlabel('Centered longitude (Japan)')
ylabel('Centered latitude (Japan)')
axis([-8 8 -10 10])

figure;
flag1 = Catalog.flag==1;
scatter(Catalog.xx(flag1), Catalog.yy(flag1), 5, log(sum(dGauss)), 'filled'); colorbar
title('Log(Sum of received background probability density) per event')
xlabel('Centered longitude (Japan)')
ylabel('Centered latitude (Japan)')
axis([-8 8 -10 10])

figure;
scatter(bwd, sum(dGauss,2), 5, 'blue', 'filled')
title('Distributed background probability density Vs. bandwidth')
xlabel('Event bandwidth')
ylabel('Distributed background probability density')

figure
idxFlag1 = find(Catalog.flag(notFlag2)==1);
scatter(bwd(idxFlag1), sum(dGauss)', 5, 'blue', 'filled')
title('Received background probability density Vs. bandwidth')
xlabel('Event bandwidth')
ylabel('Received background probability density')

sumR = sum(dGauss,2);
sumC = sum(dGauss);

%%
figure;
scatter(Catalog.xx(isFlag1), Catalog.yy(isFlag1), 5, log(lamPure), 'filled')
title('Log(Triggering event rates)')
xlabel('Centered longitude (Japan)')
ylabel('Centered latitude (Japan)')

figure;
scatter(Catalog.xx(isFlag1), Catalog.yy(isFlag1), 5, Catalog.prob(isFlag1), 'filled')
title('Background probability')
xlabel('Centered longitude (Japan)')
ylabel('Centered latitude (Japan)')

figure;
scatter(bwd(idxFlag1), Catalog.prob(isFlag1)', 5, 'blue', 'filled')
title('Background probability Vs. bandwidth')
xlabel('Event bandwidth')
ylabel('Background probability')





