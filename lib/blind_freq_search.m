function [y,deltaf,pnoise] = blind_freq_search(x,symrate,mn,M,bs,F1,F2)
%
%   Very easy to generate cycle slip
%

test_phas = (0:M-1)/M*pi/2;

corse_step = 10e6;
fine_step = 1e6;

x = DspAlg.Normalize(x,mn);
y = zeros(size(x));

tic
for pol = 1:size(x,2)
    xp = reshape(x(:,pol),bs,[]);
    
    corse_freq = -F1:corse_step:F1;
    corse_freq_phas = corse_freq*2*pi/symrate;
    
    for b = 1:size(xp,2)
        data_mat = repmat(xp(:,b),1,M);
        phas_mat = exp(1j*repmat(test_phas,bs,1));
        
        idx_vec = (1:bs) + (b-1)*bs;
        
        for k = 1:length(corse_freq)
            
            % make a frequency matrix
            phi_vec = corse_freq_phas(k) * idx_vec;
            freq_mat = exp(1j * repmat(phi_vec(:),1,M));
            
            % multiply all the matrices
            xp_tmp = data_mat.* phas_mat.* freq_mat;
            
            % search for the min distance
            xd = DspAlg.slicer(xp_tmp,mn);
            dist_phi = mean(abs(xp_tmp-xd).^2);
            [dist_corse(k) idx_phi(k)] = min(dist_phi);
        end
        [~,idx_freq] = min(dist_corse);
        rot_phas_c(b) = test_phas(idx_phi(idx_freq));
        rot_freq_c(b) = corse_freq(idx_freq);
    end
    
%     corse_freq_est = mean(rot_freq_c);
    corse_freq_est = rot_freq_c;
%     figure; plot(corse_freq_est,'.-');

    
    % run block by block
    for b = 1:size(xp,2)
        
        % fine frequency search grid
        fine_freq = corse_freq_est(b) + (-F2:fine_step:F2);
        fine_freq_phas = fine_freq*2*pi/symrate;
        
        data_mat = repmat(xp(:,b),1,M);
        phas_mat = exp(1j*repmat(test_phas,bs,1));
        
        % frequency offset index vector
        idx_vec = (1:bs);% + (b-1)*bs;
        
        for k = 1:length(fine_freq)
            
            % make a frequency matrix
            phi_vec = fine_freq_phas(k) * idx_vec;
            freq_mat = exp(1j * repmat(phi_vec(:),1,M));
            
            % multiply all the matrices
            xp_tmp = data_mat.* phas_mat.* freq_mat;
            
            % search for the min distance
            xd = slicer(xp_tmp,mn);
            dist_phi = mean(abs(xp_tmp-xd).^2);
            [dist_fine(k),idx_phi(k)] = min(dist_phi);
            
        end
        
        % process the data
        [~,idx_freq] = min(dist_fine);
        rot_phas(b) = test_phas(idx_phi(idx_freq));
        rot_freq(b) = fine_freq(idx_freq);
        tmp(:,b) = rot_freq(b) * idx_vec.'*2*pi/symrate;
    end
    
    figure; plot(tmp(:),'.-');
    
    pn = ones(bs,1) * unwrap(rot_phas*4)/4;
    df = ones(bs,1) * rot_freq;
    y(:,pol) = xp(:).* exp(1j*pn(:)).* exp(1j*unwrap(4*tmp(:))/4);
    
    pnoise(:,pol) = pn(:);
    deltaf(:,pol) = df(:);
end
toc

return


