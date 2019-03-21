function G = GjacDiffDrive(x, u)
% GjacDiffDrive adaptation to SLAM: output the jacobian of the dynamics. Returns the G matrix
%
% syms d phi theta x y
% 
% xp = x -d/phi*(sin(theta)-sin(theta+phi));
% % or xp = x +d*(cos(theta));
% yp = y +d/phi*(cos(theta)-cos(theta+phi));
% % or yp = y +d*(sin(theta));
% thetap = theta + phi;
% 
% jacobian([xp; yp; thetap], [x;y;theta])
%
%
%   INPUTS
%       x            3-by-1 vector of pose
%       u            2-by-1 vector [d, phi]'
%
%   OUTPUTS
%       G            3x3 jacobian matrix
%
%   Cornell University
%   Autonomous Mobile Robots
%   Homework 4
%   Scher, Guy
d = u(1); phi=u(2); 
theta = x(3);
% handle the case where phi=0
G = eye(length(x));
if(phi==0)
    G(1:3,1:3) =  [ 1, 0, -d*sin(theta); ...
                    0, 1,  d*cos(theta); ...
                    0, 0,             1];
else
    G(1:3,1:3) = [ 1, 0, (d*(cos(phi + theta) - cos(theta)))/phi; ...
                   0, 1, (d*(sin(phi + theta) - sin(theta)))/phi; ...
                   0, 0,                                       1];
end

