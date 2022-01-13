function Jac = jac_dista(x)
% Jac is of size (2⋅(# of points) × # of edges)
% but here it is built as Jac.' 
%
% in this case i made one point the origin,
% so i only had to keep track of 3 points.
% remember x is of size np⋅2,
% because x-coord is in the first half and y-coord in the second half.
Jac= zeros(6,6);

% row for edge # 1 (1-2)
Jac(1,1) = 2*(x(1)-x(2));
Jac(1,2) = 2*(x(2)-x(1));
Jac(1,4) = 2*(x(4)-x(5));
Jac(1,5) = 2*(x(5)-x(4));

% row for edge # 2 (2-3)
Jac(2,2) = 2*(x(2)-x(3));
Jac(2,3) = 2*(x(3)-x(2));
Jac(2,5) = 2*(x(5)-x(6));
Jac(2,6) = 2*(x(6)-x(5));

% row for edge # 3 (1-3)
Jac(3,1) = 2*(x(1)-x(3));
Jac(3,3) = 2*(x(3)-x(1));
Jac(3,4) = 2*(x(4)-x(6));
Jac(3,6) = 2*(x(6)-x(4));

% row for edge # 4 (origin-1)
Jac(4,1) = 2*x(1);
Jac(4,4) = 2*x(4);

% row for edge # 5 (origin-2)
Jac(5,2) = 2*x(2);
Jac(5,5) = 2*x(5);

% row for edge # 6 (origin-3)
Jac(6,3) = 2*x(3);
Jac(6,6) = 2*x(6);

% Jac(1,1) = (1/sqrt((x(1)-x(2))^2 + (x(4)-x(5))^2))*2*(x(1)-x(2));
% Jac(1,2) = (1/sqrt((x(1)-x(2))^2 + (x(4)-x(5))^2))*2*(x(2)-x(1));
% Jac(1,4) = (1/sqrt((x(1)-x(2))^2 + (x(4)-x(5))^2))*2*(x(4)-x(5));
% Jac(1,5) = (1/sqrt((x(1)-x(2))^2 + (x(4)-x(5))^2))*2*(x(5)-x(4));

% Jac(2,2) = (1/sqrt((x(2)-x(3))^2 + (x(5)-x(6))^2))*2*(x(2)-x(3));
% Jac(2,3) = (1/sqrt((x(2)-x(3))^2 + (x(5)-x(6))^2))*2*(x(3)-x(2));
% Jac(2,5) = (1/sqrt((x(2)-x(3))^2 + (x(5)-x(6))^2))*2*(x(5)-x(6));
% Jac(2,6) = (1/sqrt((x(2)-x(3))^2 + (x(5)-x(6))^2))*2*(x(6)-x(5));

% Jac(3,1) = (1/sqrt((x(1)-x(3))^2 + (x(4)-x(6))^2))*2*(x(1)-x(3));
% Jac(3,3) = (1/sqrt((x(1)-x(3))^2 + (x(4)-x(6))^2))*2*(x(3)-x(1));
% Jac(3,4) = (1/sqrt((x(1)-x(3))^2 + (x(4)-x(6))^2))*2*(x(4)-x(6));
% Jac(3,6) = (1/sqrt((x(1)-x(3))^2 + (x(4)-x(6))^2))*2*(x(6)-x(4));

% Jac(4,1) = (1/sqrt((x(1)-x(4))^2))*2*x(1);
% Jac(4,4) = (1/sqrt((x(1)-x(4))^2))*2*x(4);

% Jac(5,2) = (1/sqrt((x(2)-x(5))^2))*2*x(2);
% Jac(5,5) = (1/sqrt((x(2)-x(5))^2))*2*x(5);

% Jac(6,3) = (1/sqrt((x(3)-x(6))^2))*2*x(3);
% Jac(6,6) = (1/sqrt((x(3)-x(6))^2))*2*x(6);

Jac = Jac.';

end
