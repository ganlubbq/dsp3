function y = ted_AVN(px, mode)
N = length(px);
k = 1 : N;
ex = exp(-1i .* (k - 1) .* pi ./ 2);

if strcmpi(mode, 'AVN')
    s = sum(abs(px) .* ex.');
elseif strcmpi(mode, 'SLN')
    s = sum(abs(px).^2 .* ex.');
elseif strcmpi(mode, 'LOGN')
    f = log(1 + 10.*abs(px).^2);
    s = sum(f .* ex.');
elseif strcmpi(mode, 'FLN')
    s = sum(abs(px).^4 .* ex.');
end
y = -angle(s) / 2 / pi;
return
