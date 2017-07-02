% Using proportional alone doesn't work well for systems with inertia, causing oscillation and so on

clear

% static target, e.g. distance, velocity, time, frequency
target = 20;

% parameter to be controlled, using the simplest signal model
signal = 10;

gain = 0.1;

for ii = 2 : 100
  err(ii-1) = target - signal(ii-1);
  signal(ii) = signal(ii-1) + gain * err(ii-1);
end

figure;
subplot(211); plot(signal); grid on; xlabel('Iteration'); ylabel('Signal');
subplot(212); plot(err); grid on; xlabel('Iteration'); ylabel('Error');
