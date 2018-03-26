maxtrials = 10000;
N_axis = 1:maxtrials;

P = zeros(maxtrials,6);
mean_number = zeros(1,maxtrials);

for N = 1:maxtrials,
    freq = throw_die_N_times(N);
    P(N,:) = freq/N;
    mean_number(N) = mean(freq);
end;

figure(1); clf;
title('Relative frequencies of numbers shown by a die');

for k = 1:6,
    subplot(1,6,k);
    plot(N_axis,P(:,k));
    hold on;
    plot([1,maxtrials],[1/6,1/6],'r:');
    xlabel('N');
    ylabel(sprintf('f(%i)',k));
end;

figure(2); clf;
title('Mean number of points shown by a die');
plot(N_axis,mean_number);
hold on:
plot([1,maxtrials]',[3.5,3.5],'r:');

for k = 1:6,
    fprintf(1,'f(%i,%i) = %8.6f\n',k,maxtrials,P(maxtrials,k));
end;
fprintf(1,'Mean: %8.6f\n',mean_number(maxtrials));

