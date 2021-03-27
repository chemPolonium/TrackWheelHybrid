% using penalty function method
% restrain function need to be zero
% using the penalty function method
% penalty coffecient times restrain functions
% then add to the equation
% this will make the optimizer to minimize the restrain functions

% restrain functions
% let the speed to be constant
% set total speed to 15 m/s
% we need to find the best slip ratio for the vehicle
% to ride in a defined circle
v = 15;
omegaz = 0.3;
% omegaz = 1e-5;
m = 3e3;
% vrfunc: 
%   v: speed
%   r: restrain
%   func: function
% vrfunc = @(x) (norm(x)^2 - v^2)^2;
% penalty coffecient
% p = 1e4;

% vx * omegaz = Fy
% vy * omegaz = -Fx

% - vx * omegaz + Fy = 0
% vy * omegaz + Fx = 0

omegal = 50;
omegar = -50;

% x = [beta omegal omegar]
% beta: slip angle of CG
% optimfunc = @(x) norm(hybridforce(...
%     v*cos(x(1)),...
%     v*sin(x(1)),...
%     omegaz,...
%     x(2),...
%     x(3)) + [-v*cos(x(1))*omegaz*m v*sin(x(1))*omegaz*m 0])^2 + norm(x);
% x0 = [0;50;50];
% [x,fval] = fminsearch(optimfunc,x0)

optimfunc = @(x) norm(hybridforce(...
    x(1),...
    x(2),...
    x(3),...
    omegal,...
    omegar) + [-x(1)*x(3)*m x(2)*x(3)*m 0])^2;
x0 = [0;0;0];
[x,fval] = fminsearch(optimfunc,x0)