%% Path finding in a grid based graph with the dijkstra's algortihm
% obstacle is realized by varying the edge weight
function tab_dist = PathFindingDijkstra(N, source, destination, obsratio, allow_diagonal_move, allow_random_obstacle)
if nargin < 1
    N = 10;
end
if nargin < 2
    source = [1, 1];
end
if nargin < 3
    destination = [N, N];
end
if nargin < 4
    obsratio = 0.2;
end
if nargin < 5
    allow_diagonal_move = 1;
end
if nargin < 6
    allow_random_obstacle = 1;
end

nodes = ones(N);
tab_c = zeros(N + 1);
map_b = [1.0  1.0  1.0;     % white 0 - unvisited
    0.5  0.5  0.5;          % gray  1 - visited
    1.0  0.0  0.0;          % red   2 - final path
    0  0  0];               % black 3 - obstacles

% table to store the status of visited nodes
tab_vstd = zeros(N);
% table to store the previous node coordinates (i, j)
tab_prev = cell(N);
% table to store the current distance
tab_dist = ones(N) * inf;

% table to store the obstacles, using color code 3 to select the black
if allow_random_obstacle
    obstndx = randperm(N^2, obsratio * N^2);
else
    obstndx = [45, 46, 47, 55, 65, 75];
end
tab_vstd(obstndx) = 3;

% drawing
colormap(map_b);
% add one more row and column
tab_c(1:N, 1:N) = tab_vstd;
h = pcolor(tab_c);
shading faceted
axis ij

% initial
tab_dist(source(1), source(2)) = 0;

% edges in grid based graph
if allow_diagonal_move
    offset = {[-1 0], [1 0], [0 -1], [0 1], [-1 -1], [1 -1], [-1 1], [1 1]};
    % diagonal move is slower than the other directions
    weight = [1, 1, 1, 1, 1.4, 1.4, 1.4, 1.4];
    num_neighbor = 8;
else
    offset = {[-1 0], [1 0], [0 -1], [0 1]};
    weight = [1, 1, 1, 1];
    num_neighbor = 4;
end

pointer = source;

tic;
while ~ isempty(pointer)
    % marke the current node visited
    tab_vstd(pointer(1), pointer(2)) = 1;
    
    % drawing
    tab_c(1:N, 1:N) = tab_vstd;
    pcolor(tab_c);
    axis ij
    
    % find the neighbor of current node
    for ii = 1 : num_neighbor
        % validate first
        neighbor = pointer + offset{ii};
        if neighbor(1) == 0 || neighbor(2) == 0 || neighbor(1) > N || neighbor(2) > N || tab_vstd(neighbor(1), neighbor(2)) == 1
            % out of map OR has been visited
        else
            tempndx = sub2ind(size(nodes), neighbor(1), neighbor(2));
            % in the map AND unvisited AND not an obstacle
            if isempty(find(obstndx == tempndx, 1))
                dist = tab_dist(pointer(1), pointer(2)) + weight(ii);
                if dist < tab_dist(neighbor(1), neighbor(2))
                    tab_dist(neighbor(1), neighbor(2)) = dist;
                    tab_prev{neighbor(1), neighbor(2)} = pointer;
                end
            end
        end
    end
    
    % find the global min distance in unvisited set
    unvndx = find(tab_vstd == 0);
    [temp, mn] = min(tab_dist(unvndx));
    [sub1, sub2] = ind2sub(size(nodes), unvndx(mn));
    
    % update the current node
    pointer = [sub1, sub2];
    
    pause(0.01);
end

% push the destination to the path and register it as the current node
path = destination;
pointer = destination;

while ~ isempty(tab_prev{pointer(1), pointer(2)})
    % append the previous node of destination to the path and then assign
    % it to be the current node
    path = [path; tab_prev{pointer(1), pointer(2)}];
    pointer = tab_prev{pointer(1), pointer(2)};
    
    % drawing, using color code 2 to select red
    intpath = sub2ind(size(nodes), path(:,1), path(:,2));
    tab_vstd(intpath) = 2;
    tab_c(1:N, 1:N) = tab_vstd;
    pcolor(tab_c);
    axis ij
    
    pause(0.01);
end
toc;
return
