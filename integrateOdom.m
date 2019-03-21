function [x_k] = integrateOdom(x_0, u)
% integrateOdom adaptation to SLAM: Given a known initial configuration (x; y; teta) of the robot within a 
% global frame, the distance traveled d and the angle the robot turned phi,
% compute the new configuration of the robot.
% 
%   [newPose] = integrateOdom(pose_i, dist, ang_phi) returns the 
%   new configuration of the robot
% 
%   INPUTS
%       x_0      initial pose [x y teta]'
%       u        [dist; phi] vector of distance travelled, vector of angle it rotated since the sensor was last read.
%       
% 
%   OUTPUTS
%       newPose     final pose [x y teta]' x N+1
% 
% 
%   Cornell University
%   Autonomous Mobile Robots
%   Homework #2
%   Scher, Guy
% from the book probablistic robotics p.101
dist = u(1, :); phi = u(2, :);
x_k = zeros(length(x_0), length(dist)+1);
x_k(:,1) = x_0;

for i=2:(length(dist)+1)
    % the real equations "deal" with no angular rate, but matlab needs
    % help with the NaNs
    if(phi(i-1) == 0)
        x_k(1,i) = x_k(1,i-1) + dist(i-1)*cos(x_k(3,i-1)); %
        x_k(2,i) = x_k(2,i-1) + dist(i-1)*sin(x_k(3,i-1)); %
    else
        VoW = dist(i-1)/phi(i-1);
        x_k(1,i) = x_k(1,i-1) - VoW*(sin(x_k(3,i-1)) - sin(x_k(3,i-1)+phi(i-1))); %
        x_k(2,i) = x_k(2,i-1) + VoW*(cos(x_k(3,i-1)) - cos(x_k(3,i-1)+phi(i-1))); 
    end
    x_k(3,i) = x_k(3,i-1) + phi(i-1);
end
x_k(4:end,end) = x_0(4:end); %the map features don't move
