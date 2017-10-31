%Maximizing real part of a rotating complex number

clear

nsample = 1000;
s = 1; % is the target
x = exp(1i * pi / 4);

phi(1) = 0; % is the initial value of parameter to be controlled
r(1) = 0;
for ii = 1 : 100
    r(ii + 1) = x * exp(-1i * phi(ii));
    d(ii) = r(ii + 1) - r(ii);
end
