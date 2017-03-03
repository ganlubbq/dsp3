%% Path finding demo with the dijkstra's algortihm

function PathFindingDijkstra(N, source, destination)
if nargin < 1
    N = 10;
end
if nargin < 2
    source = [1, 1];
end
if nargin < 3
    destination = [N, N];
end

tic;

nodes = ones(N);
tab_c = ones(N + 1);
map_c = [1.0  1.0  1.0;
    0.5  0.5  0.5];

% table to store the status of visited nodes
tab_vstd = zeros(N);
% table to store the previous node coordinates (i, j)
tab_prev = cell(N);
% table to store the current distance
tab_dist = ones(N) * inf;

% drawing
tab_c(1:N, 1:N) = tab_vstd;
h = pcolor(tab_c);
colormap(map_c);
shading faceted
axis ij

% initial
tab_dist(source(1), source(2)) = 0;
% edges
offset = {[-1 0], [1 0], [0 -1], [0 1], [-1 -1], [1 -1], [-1 1], [1 1]};
% diagonal move is slower than the other directions
weight = [1, 1, 1, 1, 1.4, 1.4, 1.4, 1.4];

pointer = source;

while ~ isempty(pointer)
    % marke the current node visited
    tab_vstd(pointer(1), pointer(2)) = 1;
    
    % drawing
    tab_c(1:N, 1:N) = tab_vstd;
    pcolor(tab_c);
    axis ij
    
    % find the neighbor of current node
    for ii = 1 : 8
        % validate first
        neighbor = pointer + offset{ii};
        if neighbor(1) == 0 || neighbor(2) == 0 || neighbor(1) > N || neighbor(2) > N
            % out of map
        else
            % in the map, make sure it is unvisited
            if tab_vstd(neighbor(1), neighbor(2)) == 0
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

path = destination;
tab_path = ones(N);

destndx = sub2ind(size(nodes), destination(1), destination(2));

while ~ isempty(tab_prev{destndx})
    path = [path; tab_prev{destndx}];
    pointer = tab_prev{destndx};
    destndx = sub2ind(size(nodes), pointer(1), pointer(2));
    
    intpath = sub2ind(size(nodes), path(:,1), path(:,2));
    tab_path(intpath) = 0;
    pcolor(~tab_path);
    axis ij
    
    pause(0.01);
end
toc;
return
