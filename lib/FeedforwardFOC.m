function y = FeedforwardFOC(params)
%FREQ_OFFSET frequency offset estimation routine
% using bell-lab style
%
% CopyRight:Wang Dawei EIE PolyU   $Date:16/3/2010
symbol_rate = params.baudrate;
x = params.in;
y = x;
if params.active
    disp('Begin FOE (frequency offset estimation)...')
    N = length(x);
    data2 = x(2:end,:).^4;
    data1 = x(1:end-1,:).^4;
    data  = data2.*conj(data1);
    delta_sum = sum(data);
    angle_acq = angle(delta_sum)/4;
    dphi = -mean(angle_acq);
    fo = (1:N) * dphi;
    y = x.*(exp(1j*fo.')*ones(1,size(x,2)));

    deltaf = dphi*symbol_rate/(2*pi)
    disp('Done')
end
end