clear

% connect_table = {[2,3,4,7],[1,3,5,8],[1,2,6,9],[1,5,6,7],[2,4,6,8],[3,4,5,9],[1,4,8,9],[2,5,7,9],[3,6,7,8]};
connect_table = {[2,3,5,8],[1,4],[1,4,6],[2,3,6,7],[1,6,8,9],[3,4,5,7,9],[4,6,9],[1,5,9],[5,6,7,8]};

% number of vertex
N = 9;

% number of color
C = 3;

% randomly initialize the color
color_table = floor(rand(N,1) * C + 1);

% initialize the cost
F = 0;
for cc = 1 : C
    ndxC = find(color_table == cc);
    
    % at least we should have one conflict
    if length(ndxC) < 2
        continue;
    else
        for ii = 1 : length(ndxC)
            adj_node = connect_table{ndxC(ii)};
            for jj = ii + 1 : length(ndxC)
                if find(adj_node == ndxC(jj), 1)
                    F = F + 1;
                end
            end
        end
    end
end

% the adjacent-color table
adjcolor_table = zeros(N, C);

% initialize the adjacent-color table
for ii = 1 : N
    adj_node = connect_table{ii};
    colorndx = color_table(adj_node);
    for jj = 1 : C
        adjcolor_table(ii, jj) = sum(colorndx == jj);
    end
end

% initialize tabu table
tabu_table = zeros(N, C);

% tabu length
tt = 2 * F;

% iteration counter
ic = 0;

% main iteration
while F > 0
    
    % minus all positive tabu by 1
    tabu_table(tabu_table > 0) = tabu_table(tabu_table > 0) - 1;
    
    % find the best one-move
    for ii = 1 : N
        conflict(ii) = adjcolor_table(ii, color_table(ii));
    end
    
    ndxConf = find(conflict > 0);
    
    ndxN = ndxConf(floor(rand() * length(ndxConf) + 1));
    
    ndxi = color_table(ndxN);
    
    for ndxj=1 : C
        % if there is no tabu, make the move
        if (tabu_table(ndxN,ndxj) == 0) && (ndxj ~= ndxi)
            
            color_table(ndxN) = ndxj;
            
            tabu_table(ndxN, ndxi) = tt;
            
            maxdelta = adjcolor_table(ndxN, ndxi) - adjcolor_table(ndxN, ndxj);
            F = F - maxdelta;
            
            % update the adjacent-color table
            adjNode = connect_table{ndxN};
            adjcolor_table(adjNode, ndxi) = adjcolor_table(adjNode, ndxi) - 1;
            adjcolor_table(adjNode, ndxj) = adjcolor_table(adjNode, ndxj) + 1;
            break;
        end
    end
    ic = ic + 1;
end
ic
open color_table
