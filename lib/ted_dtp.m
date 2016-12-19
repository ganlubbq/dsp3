function y = ted_dtp(px)
%Delay tap timing error detector

x = abs(px(1:2:end-1)).^2;
y = abs(px(2:2:end)).^2;

% cosine
distx = (x-y)./sqrt(2);

% sine
% distx = (y-x)./sqrt(2);

y = mean(distx);


% rx = real(px(1:2:end-2));
% ry = real(px(3:2:end));
% ix = real(px(1:2:end-2));
% iy = real(px(3:2:end));
% % distx = (x-y)./sqrt(2);
% distr = (abs(ry)+abs(rx)).^2./sqrt(2);
% disti = (abs(iy)+abs(ix)).^2./sqrt(2);
% y = mean(distr+disti);