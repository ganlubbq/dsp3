function b = sym2bit(x,mn)
%SYM2BIT Convert gray mapping symbols to binary bits

map16 = [...
    0 0 0 0; 0 0 0 1; 0 0 1 1; 0 0 1 0;
    0 1 0 0; 0 1 0 1; 0 1 1 1; 0 1 1 0;
    1 1 0 0; 1 1 0 1; 1 1 1 1; 1 1 1 0;
    1 0 0 0; 1 0 0 1; 1 0 1 1; 1 0 1 0 ];
map4 = [...
    0 0; 0 1;
    1 0; 1 1 ];

switch mn
    case 2
        b = x*2 - 1;
    case 4
        b = map4(x,:);
    case 16
        b = map16(x,:);
end