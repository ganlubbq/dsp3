function CheckEvenNargin( varargin )
%CHECKEVENNARGIN Summary of this function goes here
%   Detailed explanation goes here

moda = mod( length(varargin), 2 );

if moda
    error('CHECKEVENNARGIN::Number of inputs should be EVEN')
end

return

