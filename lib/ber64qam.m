function BER = ber64qam(xin,ref)
%BER16QAM_EXP Bit error rate test for 16-QAM experiment
%   diff: flag of differential decoding
%   Vectoring method
%
%   Copyright: WANG Dawei [EIE PolyU]   $Date:16/3/2010$

% Decision
% data = decision_qam3(xin(1000:end),6).';

% load(filename);

local_max = Hist2(xin);
ref_const = local_max(:);

% figure(99); 
% plot(xin,'.');          hold on
% plot(ref_const,'ro');   hold off

map = decision_qam3(ref_const,6);

data = zeros(size(xin));
for kk = 1:length(xin)
    [~,idx] = min(abs(xin(kk)-ref_const));
    data(kk) = map(idx);
end

data = data.';


reflen = length(ref(:,1));
a = mod(length(data),reflen);
data = data(a:end,:);
N   = length(data)-1;

L = 1000;
ref1_tmp = ref(1:L,1);
data_Test_bit = differential_decode(data);
[c1,lags1] = xcorr(ref1_tmp,data_Test_bit(1:reflen,1));
[~,idx] = max(c1);

dataTest_bit = circshift(data_Test_bit,lags1(idx));

ref_Data = repmat(ref,N/reflen,1);
BER =nnz(dataTest_bit-ref_Data)/N/6;

return

