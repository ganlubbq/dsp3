%% Demo program to draw a shaded and squared chessboard in matlab
% one can play chesses by varying the color matrix elements

clear;

N = 100;

% empty chessboard with all white
board = ones(N+1);

% play
tic
for rr = 1:N
    board(rr,:) = rand(1,N+1) > 0.5;
    
    % draw
    h = pcolor(board);
    
    colormap(gray(2));
    
    shading faceted
    
    axis ij
    axis square
    
    pause(0.01)
end
toc
