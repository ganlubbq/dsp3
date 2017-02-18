% Bit error rate test using anti-rotation method
% Possible reference: x, y, -x, -y
%
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Default

function [ser, xout, pos] = ber_anti_rotation(xin, ref_sym)

% Init, ref1 for x and ref2 for y
ref1 = ref_sym(:,1);
ref2 = ref_sym(:,2);

% -x and -y
ref1_inv = -ref1;
ref2_inv = -ref2;

% data length should be multiple integer of reference length
ref_len = length(ref1);
data_len = length(xin);
data = [real(xin),imag(xin)];
% if not, discard few data in the end
data = data(1:end-mod(data_len,ref_len),:);
data_len = length(data(:,1));

% Multi-level Decision for 16-qam
thr1 = -2;
thr2 = 0;
thr3 = 2;
data(data>=thr3) = 3;
data(data>=thr2 & data<thr3) = 1;
data(data>=thr1 & data<thr2) = -1;
data(data<thr1) = -3;

ref1_tmp = ref1(1:1000);
ref2_tmp = ref2(1:1000);
ref1_inv_tmp = ref1_inv(1:1000);
ref2_inv_tmp = ref2_inv(1:1000);

[c1,lags1] = xcorr(ref1_tmp,data(1:ref_len,1));
[c2,lags2] = xcorr(ref2_tmp,data(1:ref_len,1));
[c3,lags3] = xcorr(ref1_inv_tmp,data(1:ref_len,1));
[c4,lags4] = xcorr(ref2_inv_tmp,data(1:ref_len,1));
[~,idx] = max([max(c1),max(c2),max(c3),max(c4)]);

switch idx
    case 1
        refi = ref1;
        [~,index] = max(c1);
        data_tmp = circshift(data,lags1(index));
        lags1(index)
        c1 = xcorr(ref2_tmp,data_tmp(1:ref_len,2));
        c2 = xcorr(ref2_inv_tmp,data_tmp(1:ref_len,2));
        if  max(c1) > max(c2)
            refq = ref2; 
            pos = 1;
            xout = data_tmp;
        else
            refq = ref2_inv;
            pos = 2;
            %xout = data_tmp;%
            xout = [data_tmp(:,1),-data_tmp(:,2)];
        end
        refi = repmat(refi,data_len/ref_len,1);
        refq = repmat(refq,data_len/ref_len,1);
        symerr_num_i = nnz(refi - data_tmp(:,1));
        symerr_num_q = nnz(refq - data_tmp(:,2));
        
    case 2
        refi = ref2;
        [~,index] = max(c2);
        data_tmp = circshift(data,lags2(index));
        lags2(index)
        c1 = xcorr(ref1_tmp,data_tmp(1:ref_len,2));
        c2 = xcorr(ref1_inv_tmp,data_tmp(1:ref_len,2));
        if  max(c1) > max(c2)
            refq = ref1;
            pos = 3;
            %xout = data_tmp;%
            xout = [data_tmp(:,2),data_tmp(:,1)];
        else
            refq = ref1_inv;
            pos = 4;
            xout = [-data_tmp(:,2),data_tmp(:,1)];
        end
        refi = repmat(refi,data_len/ref_len,1);
        refq = repmat(refq,data_len/ref_len,1);
        symerr_num_i = nnz(refi - data_tmp(:,1));
        symerr_num_q = nnz(refq - data_tmp(:,2));
        
    case 3
        refi = -ref1;
        [~,index] = max(c3);
        data_tmp = circshift(data,lags3(index));
        lags3(index)
        c1 = xcorr(ref2_tmp,data_tmp(1:ref_len,2));
        c2 = xcorr(ref2_inv_tmp,data_tmp(1:ref_len,2));
        if  max(c1) > max(c2)
            refq = ref2;
            pos = 5;
            %xout = data_tmp;%
            xout = [-data_tmp(:,1),data_tmp(:,2)];
        else
            refq = ref2_inv;
            pos = 6;
            xout = [-data_tmp(:,1),-data_tmp(:,2)];
        end
        refi = repmat(refi,data_len/ref_len,1);
        refq = repmat(refq,data_len/ref_len,1);
        symerr_num_i = nnz(refi - data_tmp(:,1));
        symerr_num_q = nnz(refq - data_tmp(:,2));
        
    case 4
        refi = -ref2;
        [~,index] = max(c4);
        data_tmp = circshift(data,lags4(index));
        lags4(index)
        c1 = xcorr(ref1_tmp,data_tmp(1:ref_len,2));
        c2 = xcorr(ref1_inv_tmp,data_tmp(1:ref_len,2));
        if  max(c1) > max(c2)
            refq = ref1;
            pos = 7;
            xout = [data_tmp(:,2),-data_tmp(:,1)];
        else
            refq = ref1_inv;
            pos = 8;%
            %xout = data_tmp;
            xout = [-data_tmp(:,2),-data_tmp(:,1)];
        end
        refi = repmat(refi,data_len/ref_len,1);
        refq = repmat(refq,data_len/ref_len,1);
        symerr_num_i = nnz(refi - data_tmp(:,1));
        symerr_num_q = nnz(refq - data_tmp(:,2));
end

ser = (symerr_num_i + symerr_num_q) / data_len;

return

