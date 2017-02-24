%% TEST SCRIPT FOR ESTIMATING CARRIER PHASE USING MLSE
% IMPLEMENTING MLSE ESTIMATION BY VITERBI ALGORITHM
%
% complied 2015-12-19
%
%% TODO: add periodic pilot

clear

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1234));

format compact

bitpersym = 2;

% number of states in this case
M = 2^bitpersym;

%------------global setting-----------
% memory depth, number of consecutive symbols for metric calculation...
% v = 4;

% truncated survivor length...
% decision delay, any stage before T-L is forced to make decision... L
% should be at least 5 times larger than memory depth to ensure paths merge
% before T-L...
L = 20;

% segment length
segLen = 256;

% total symbol length
symlen = 2^15 +L;

% differential decoding switch
dd = 1;

% loop delay
delay = 1;

% combined laser linewidth
lsrlw = 1e6;

% symbol rate
symrate = 28e9;

% osnr
osnr = 11:16;

% constellation
C = constellation(M);

% mean power
P = sum(abs(C).^2) / M;

%
delta = lsrlw / symrate;

%------------data path--------------
%
fprintf('generating binary data...\n');
bitTx = randi([0 1], bitpersym, symlen);

%
fprintf('symbolizing...\n');
symTx = symbolizerGrayQam(bitTx);

% laser noise
lsrPNstd = sqrt(2*pi*lsrlw/symrate);
lsrPNrnd = rand(1,symlen);
lsrPN = 0 + cumsum((lsrPNrnd-mean(lsrPNrnd))/calcrms(lsrPNrnd)*lsrPNstd);

%
Tx = symTx.' .* exp(1i * lsrPN);
Ts = symTx(1:L).';

%----------start testing----------
%
for iOSNR = 1:length(osnr)
    
    fprintf('starting osnr loop #%d...\n', iOSNR);
    
    % add noise
    % Rx = setosnr(Tx,osnr(iOSNR),symrate);
    esno = 10^(osnr(iOSNR)*0.1)*(12.5e9)/symrate;
    pn = sum(abs(Tx).^2)/length(Tx)/esno;
    noise = genWGN(size(Tx,1),size(Tx,2),pn,'linear','complex');
    Rx = Tx+noise;
    
    % truncate survivor path (pointer list) with length L
    % only M survivors is kept in mind
    Srv = zeros(M,L);
    
    %0   1   2   3   4   5   6   7   8   9   10   11
    %-----------------------------------------------
    %|   processed data  |      M survivors        |
    %-----------------------------------------------
    % time 11 is the current stage
    
    % distance matrix (score)
    % each of the M states in current step has M predecessors
    % first M stands for current stage, second M stands for previous stage
    U = zeros(symlen,M,M);
    
    % decisions
    Dec = zeros(symlen,1);
    
    % phase noise estimations
    estPN = zeros(symlen,1);
    
    % training sequence to start up (solving the phase ambiguity)
    % the value stored in time i is the survivor index in i-1
    for iTmp = 2:L
        Srv(:,iTmp) = find(Ts(iTmp-1)==C);
    end
    
    %---------process initial symbols----------
    % time index of overall
    kk = L+1;
    
    % time index of useful and current stage
    ndx = 1;
    
    % symbol index of last training symbol
    ii = find(Ts(L)==C);
    
    % using past L symbols to estimate phase noise
    metric = sum(conj(Rx(kk-L:kk-1)).*Ts)/(P*L);
    estPN(1) = angle(metric);
    
    % make first decision before shift srv (moving forward)...
    Dec(1) = Ts(1);
    Srv = circshift(Srv,[0 -1]);
    
    % initializing U & updating srv
    for jj = 1:M
        
        % Srv(:,L), last column of srv, stores the srv pointer for current stage
        % for initializing purpose, this pointer points to the last symbol of Ts
        Srv(jj,L) = ii;
        
        % calculating the distance between the previous stage and current
        % stage...at initialization, the previous srv is the last symbol of Ts
        % for all states in current stage
        U(ndx,jj,ii) = -abs(Rx(kk).*metric-C(jj))^2;
    end
    
    % flush memory, always keep L recent symbols as memory for M states...as
    % can be shown below, memory can be obtained by backtracking the survivor
    % path by depth of L...
    % TODO: L should be the memory length for memory channel...
    memory = zeros(M,L);
    
    % flush metric
    % branch metric for current stage...
    metric = zeros(M,1);
    
    %--------process rest of symbols----------
    %
    fprintf('processing...\n');
    for kk = L+2:length(Rx)
        
        % time index of useful and current stage
        ndx = kk - L;
        
        for ii = 1:M
            iRow = ii;
            
            % backtrack srv to update memory for each state
            for mm = L-1:-1:1
                memory(ii,mm) = C(Srv(iRow,mm+1));
                iRow = Srv(iRow,mm+1);
            end
            
            % M possible paths up to time i-1 (previous stage)
            memory(ii,L) = C(ii);
            
            % using past L symbols as memory to estimate phase noise
            metric(ii) = sum(conj(Rx(kk-L:kk-1)) .* memory(ii,:))/(P*L);
        end
        
        % make decision before shift srv (move forward)...by truncating the
        % survivor length, everything before survivor will be forced as
        % firm decision...before moving forward, the oldest survivor will become
        % a decision...
        if length(unique(Srv(:,2))) == 1
            % if all states are pointing to the same previous state,
            % that is the decision...
            Dec(ndx) = C(unique(Srv(:,2)));
        else
            % in order to address the error event (path diverge), multiple rules
            % apply such as random pick, max metric, majority vote...
            Dec(ndx) = C(Srv(1,2));
            % TODO: add periodic pilot symbols, so that paths will converge at
            % pilots
        end
        
        Srv = circshift(Srv, [0 -1]);
        
        % updating U & srv according to Viterbi algorithm
        % jj stands for current stage, total MxM combinations
        for jj = 1:M
            for ii = 1:M
                U(ndx,jj,ii) = U(ndx-1,ii,Srv(ii,L-1)) - abs(Rx(kk).*metric(ii)-C(jj)).^2;
            end
            % find the survivor in previous stage
            [maxU,ndxMU] = max(U(ndx,jj,:));
            Srv(jj,L) = ndxMU;
            estPN(ndx) = angle(metric(ndxMU));
        end
    end
    
    % figure;
    % plot(lsrPN); hold;
    % plot(-estPN,'g'); hold off;
    
    bitRx = slicerGrayQam(Dec,M);
    
    bec(iOSNR) = sum(sum(bitTx(:,1:end-L) ~= bitRx(:,1:end-L)));
    ber(iOSNR) = bec(iOSNR) / symlen / bitpersym;
    
    fprintf('process returns with ber %.2e\n', ber(iOSNR));
end

opt = T_BER_mQAM(osnr,M,symrate);

figure;
plot(osnr,log10(opt)); hold;
plot(osnr,log10(ber),'s-'); grid on; hold off;

