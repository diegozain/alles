function A = adj_grid(n,m)

%					z
%			-----------------
%
%		|	1-------a-------3
%		|	| o-----------o |
%		|	| |           | |
%	x	|	b |           | d
%		|	| |           | |
%		|	| o-----------o |
%		|	2-------c-------4
%
%
%   grid inside -o- is inner grid.
%   -o- is inner grid boundary.

%----------------------------------------------------------------------
%             build A 
%----------------------------------------------------------------------

I = zeros(n*m,1);
J = zeros(n*m,1);
V = zeros(n*m,1);

%--------------------------------------
% inner
%--------------------------------------

k=1;
for l=2:n-1
    for h=2:m-1
        i = coordP([l,h],n);
        % vertical (down) neighbors
        j = i + 1;
        J(k) = j;
        I(k) = i;
        V(k) = 1;
        k=k+1;
        % vertical (up) neighbors
        j = i - 1;
        J(k) = j;
        I(k) = i;
        V(k) = 1;
        k=k+1;
        % horizontal (right) neighbors
        j = i + n;
        J(k) = j;
        I(k) = i;
        V(k) = 1;
        k=k+1;
        % horizontal (left) neighbors
        j = i - n;
        J(k) = j;
        I(k) = i;
        V(k) = 1;
        k=k+1;
        % % inner nodes
        % J(k) = i;
        % I(k) = i;
        % V(k) = 1;
        % k=k+1;
    end
end

%--------------------------------------
% neu edge b
%--------------------------------------

h=1;
for l=2:n-1
    i = coordP([l,h],n);
    % vertical (down) neighbors
    j = i + 1;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % vertical (up) neighbors
    j = i - 1;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % horizontal (right) neighbors
    j = i + n;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % % inner nodes
    % J(k) = i;
    % I(k) = i;
    % V(k) = 1;
    % k=k+1;
end

%--------------------------------------
% corner 1
%--------------------------------------

i = 1;
% vertical (down) neighbor
j = i + 1;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% horizontal (right) neighbor
j = i + n;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% % inner node
% J(k) = i;
% I(k) = i;
% V(k) = 1;
% k=k+1;

%--------------------------------------
% corner 2
%--------------------------------------

i = n;
% vertical (up) neighbor
j = i - 1;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% horizontal (right) neighbor
j = i + n;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% % inner node
% J(k) = i;
% I(k) = i;
% V(k) = 1;
% k=k+1;

%--------------------------------------
% robin edge a
%--------------------------------------

l=1;
for h=2:m-1
    i = coordP([l,h],n);
    % vertical (down) neighbors
    j = i + 1;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % horizontal (left) neighbors
    j = i - n;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % horizontal (right) neighbors
    j = i + n;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % % inner nodes
    % J(k) = i;
    % I(k) = i;
    % V(k) = 1;
    % k=k+1;
end

%--------------------------------------
% robin edge d
%--------------------------------------

h=m;
for l=2:n-1
    i = coordP([l,h],n);
    % vertical (down) neighbors
    j = i + 1;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % vertical (up) neighbors
    j = i - 1;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % horizontal (left) neighbors
    j = i - n;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % % inner nodes
    % J(k) = i;
    % I(k) = i;
    % V(k) = 1;
    % k=k+1;
end

%--------------------------------------
% robin edge c
%--------------------------------------

l=n;
for h=2:m-1
    i = coordP([l,h],n);
    % horizontal (left) neighbors
    j = i - n;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % vertical (up) neighbors
    j = i - 1;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % horizontal (right) neighbors
    j = i + n;
    J(k) = j;
    I(k) = i;
    V(k) = 1;
    k=k+1;
    % % inner nodes
    % J(k) = i;
    % I(k) = i;
    % V(k) = 1;
    % k=k+1;
end

%--------------------------------------
% corner 3
%--------------------------------------

i = n*(m-1) + 1;
% vertical (down) neighbor
j = i + 1;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% horizontal (left) neighbor
j = i - n;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% % inner node
% J(k) = i;
% I(k) = i;
% V(k) = 1;
% k=k+1;

%--------------------------------------
% corner 4
%--------------------------------------

i = n*m;
% vertical (up) neighbor
j = i - 1;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% horizontal (left) neighbor
j = i - n;
J(k) = j;
I(k) = i;
V(k) = 1;
k=k+1;
% % inner node
% J(k) = i;
% I(k) = i;
% V(k) = 1;

%-------------------------------------------------------------------------%

A = sparse(I,J,V);

end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

function u = coordP(XX,n)
% given grid coordinate, find number coordinate of point.
%
% XX <-- is a vector whose entries are the coordinate of the point in the 
%        grid.
% n  <-- horizontal grid size.
%
% u  <-- number coordinate of point.

	i = XX(1);
	j = XX(2);
	
	u = i + (j-1)*n;
end