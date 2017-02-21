function c = constellation(mn)
%CONSTELLATION Give the canonical form of uncoded constellations for mn-QAM
% The integer index is from topleft to bottomright by columns


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

return


