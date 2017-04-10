function c = constellation(mn)
%CONSTELLATION Give the canonical form of uncoded constellations for mn-QAM
% The integer index is from topleft to bottomright by columns

% % reference
% h = modem.qammod('M',mn);
% cr = h.Constellation;

switch mn
    case 2
        c = [-1; 1];
        
    case 4
        c = [-1+1j,-1-1j,1+1j,1-1j;];
        
    case 8
        c = [-2, -2j, 2j, 2, 3+3j, 3-3j, -3-3j, -3+3j];
        
    case 16
        ci = [-3 -3 -3 -3 -1 -1 -1 -1 1 1 1 1 3 3 3 3];
        cq = [3 1 -1 -3 3 1 -1 -3 3 1 -1 -3 3 1 -1 -3];
        c = ci + 1j*cq;
        
    case 32
        c = [-3 + 5i,-1 + 5i,-1 - 5i,-3 - 5i,-5 + 3i,-5 + 1i,-5 - 1i,-5 - 3i,-3 + 3i,-3 + 1i,-3 - 1i,-3 - 3i,-1 + 3i,-1 + 1i,-1 - 1i,-1 - 3i,1 + 3i,1 + 1i,1 - 1i,1 - 3i,3 + 3i,3 + 1i,3 - 1i,3 - 3i,5 + 3i,5 + 1i,5 - 1i,5 - 3i,3 + 5i,1 + 5i,1 - 5i,3 - 5i];
    
    case 64
        ci = [-7*ones(1,8) -5*ones(1,8) -3*ones(1,8) -1*ones(1,8) ones(1,8) 3*ones(1,8) 5*ones(1,8) 7*ones(1,8)];
        temp = [7 5 3 1 -1 -3 -5 -7];
        cq = [temp temp temp temp temp temp temp temp];
        c = ci + 1j*cq;
        
    case 256
        ci = [-15*ones(1,16) -13*ones(1,16) -11*ones(1,16) -9*ones(1,16) -7*ones(1,16) -5*ones(1,16) -3*ones(1,16) -1*ones(1,16) 1*ones(1,16) 3*ones(1,16) 5*ones(1,16) 7*ones(1,16) 9*ones(1,16) 11*ones(1,16) 13*ones(1,16) 15*ones(1,16)];
        temp = 15:-2:-15;
        cq = [temp temp temp temp temp temp temp temp temp temp temp temp temp temp temp temp];
        c = ci + 1j*cq;
        
    otherwise
        warning('unsupported modulation format'); keyboard;
end

c = c(:);

return