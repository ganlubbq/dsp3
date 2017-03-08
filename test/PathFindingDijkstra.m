%% Path finding in a grid based graph with the dijkstra's algortihm
% obstacle is realized by varying the edge weight
%
% Test Cases:
%   PathFindingDijkstra(10, [1 1], [10 10], 0, 0, 0, 0);
%   PathFindingDijkstra(10, [1 1], [10 10], 0, 0, 0, 1);
%   PathFindingDijkstra(10, [1 1], [10 10], 0, 1, 0, 0);
%   PathFindingDijkstra(10, [1 1], [10 10], 0, 1, 0, 1);
%   PathFindingDijkstra(10, [1 1], [10 10], 0.1, 0, 0, 0);
%   PathFindingDijkstra(10, [1 1], [10 10], 0.1, 0, 0, 1);
%   PathFindingDijkstra(10, [1 1], [10 10], 0.1, 1, 0, 0);
%   PathFindingDijkstra(10, [1 1], [10 10], 0.1, 1, 0, 1);
%   while 1; PathFindingDijkstra(10, [1 1], [10 10], .25, 0, 1, 0); pause(0.5); end
%   while 1; PathFindingDijkstra(10, [1 1], [10 10], .25, 0, 1, 1); pause(0.5); end
%   while 1; PathFindingDijkstra(10, [1 1], [10 10], .25, 1, 1, 0); pause(0.5); end
%   while 1; PathFindingDijkstra(10, [1 1], [10 10], .25, 1, 1, 1); pause(0.5); end
%
function [path, tab_dist] = PathFindingDijkstra(N, source, destination, obsratio, allow_diagonal_move, allow_random_obstacle, allow_heuristic)
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
if nargin < 7
    allow_heuristic = 0;
end

nodes = ones(N);
tab_c = zeros(N + 1);

map_3 = [1.0  1.0  1.0;     % white 0 - unvisited
    1.0  0.0  0.0;          % red   1 - final path
    0.5  0.5  0.5];         % gray  2 - visited

map_4 = [1.0  1.0  1.0;     % white 0 - unvisited
    0.5  0.5  0.5;          % gray  1 - visited
    1.0  0.0  0.0;          % red   2 - final path
    0  0  0];               % black 3 - obstacles

if obsratio == 0
    color_code_unvst = 0;
    color_code_vsted = 2;
    color_code_path  = 1;
    color_code_obst  = 3;
    map_c = map_3;
else
    color_code_unvst = 0;
    color_code_vsted = 1;
    color_code_path  = 2;
    color_code_obst  = 3;
    map_c = map_4;
end

% table to store the status of visited nodes
tab_vstd = ones(N) * color_code_unvst;
% table to store the previous node coordinates (i, j)
tab_prev = cell(N);
% table to store the current distance defined by edge
tab_dist = ones(N) * inf;
% table to store the current score defined by distance + heuristic
tab_scre = ones(N) * inf;

% table to store the obstacles, using color code 3 to select the black;
% force the source and goal nodes not being obstacles
if allow_random_obstacle
    obstndx = randperm(N^2, round(obsratio * N^2));
else
    if obsratio == 0
        obstndx = [];
    else
%         obstndx = [38, 48, 58, 68, 67, 66, 65, 64];
        obstndx = [36, 45, 47, 56, 57] - 1;
    end
end
tab_vstd(obstndx) = color_code_obst;
tab_vstd(source(1), source(2)) = color_code_unvst;
tab_vstd(destination(1), destination(2)) = color_code_unvst;

% drawing
colormap(map_c);
% add one more row and column
tab_c(1:N, 1:N) = tab_vstd;
h = pcolor(tab_c);
shading faceted
axis ij

% initial
tab_dist(source(1), source(2)) = 0;
tab_scre(source(1), source(2)) = 0;

% edges in grid based graph, if diagonal move is allowed, the wavefront is
% round, or curved. If diagonal move is not allowed, the wavefront is flat
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
iteration = 0;

tic;
while ~ isempty(pointer)
    % marke the current node visited
    tab_vstd(pointer(1), pointer(2)) = color_code_vsted;
    iteration = iteration + 1;
    
    % drawing
    tab_c(1:N, 1:N) = tab_vstd;
    pcolor(tab_c);
    axis ij
    title(sprintf('iteration %d', iteration));
    
    % find the neighbor of current node
    for ii = 1 : num_neighbor
        % validate first
        neighbor = pointer + offset{ii};
        if neighbor(1) == 0 || neighbor(2) == 0 || neighbor(1) > N || neighbor(2) > N || tab_vstd(neighbor(1), neighbor(2)) == color_code_vsted
            % out of map OR has been visited, skip this neighbor
        elseif tab_vstd(neighbor(1), neighbor(2)) == color_code_unvst
            % in the map AND unvisited AND not an obstacle, update its
            % score accordingly
            dist = tab_dist(pointer(1), pointer(2)) + weight(ii);
            if dist < tab_dist(neighbor(1), neighbor(2))
                tab_dist(neighbor(1), neighbor(2)) = dist;
                tab_prev{neighbor(1), neighbor(2)} = pointer;
                if allow_heuristic
                    tab_scre(neighbor(1), neighbor(2)) = dist + heuristic(destination, neighbor);
                else
                    tab_scre(neighbor(1), neighbor(2)) = dist;
                end
            end
        end
    end
    
    % find the global min scre in unvisited set; in case multiple nodes
    % have the same score, select one with least distance
    unvndx = find(tab_vstd == color_code_unvst);
    list_scre = tab_scre(unvndx);
%     [temp, mn] = min(tab_scre(unvndx));
    unvndx_minscre = unvndx(list_scre == min(list_scre));
    list_dist = tab_dist(unvndx_minscre);
    [temp, mn] = min(list_dist);
    
    % some places you simply can not go
    if isinf(temp)
        break;
    else
        [sub1, sub2] = ind2sub(size(nodes), unvndx_minscre(mn));
    end
    % update the current node
    pointer = [sub1, sub2];
    
    % if the destination is reached, terminate, otherwise the while loop
    % will run until all nodes are visited, i.e. no solution is found;
    % comment it out if you want to loop all the nodes.
    if sub1 == destination(1) && sub2 == destination(2)
        break;
    end
    
    pause(0.01);
end

% push the destination to the path and register it as the current node
path = destination;
pointer = destination;
final_score = 0;
while ~ isempty(tab_prev{pointer(1), pointer(2)})
    % append the previous node of destination to the path and then assign
    % it to be the current node
    final_score = final_score + manhattan(pointer, tab_prev{pointer(1), pointer(2)}); 
    path = [path; tab_prev{pointer(1), pointer(2)}];
    pointer = tab_prev{pointer(1), pointer(2)};
    
    % drawing, using color code 2 to select red
    intpath = sub2ind(size(nodes), path(:,1), path(:,2));
    tab_vstd(intpath) = color_code_path;
    tab_c(1:N, 1:N) = tab_vstd;
    pcolor(tab_c);
    axis ij
    title(sprintf('score of %.1f at iteration %d', final_score, iteration));
    
    pause(0.01);
end
toc;
return

% educated guess of cost from node "neighbor" to node "destination"
% it is highly recommended to UNDERESTIMATE the heuristic (coeff <= 1)
function h = heuristic(destination, neighbor)
h2 = (destination(1) - neighbor(1)) .^ 2 + (destination(2) - neighbor(2)) .^ 2;
h = 1.4 * sqrt(h2);
return

% get manhattan distance of two points
function h = manhattan(destination, neighbor)
h = abs(destination(1) - neighbor(1)) + abs(destination(2) - neighbor(2));
return
