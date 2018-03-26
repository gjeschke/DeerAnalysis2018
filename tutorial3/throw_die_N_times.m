function freq = throw_die_N_times(N)

freq = zeros(1,6);
for k = 1:N,
    points = die;
    freq(points) = freq(points) + 1;
end;