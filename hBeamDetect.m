function meas = hBeamDetect(mu)
%   hBeamDetect(mu) returns the expected measurements of the 
%                   beam sensor from each obstacle. See hw6writeup for
%                   explanations
%   INPUTS
%       mu      [x y teta x1 y1 x2 y2 ... xi yi]'
%       
%   OUTPUTS
%       meas    [measx1 measy1 .... measxi measyi]'
% 
%   Cornell University
%   Autonomous Mobile Robots
%   Homework #6
%   Scher, Guy

Nm = (length(mu) - 3)/2;
robot = mu(1:2); theta = mu(3);
T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
meas = zeros(Nm*2, 1);

for i=1:2:(Nm*2)
    feat_pos = mu( (3+i) : (3+i+1) );
    meas(i:i+1) = T*(feat_pos - robot);
end
% if it's behind me, it's definitely not in view
% meas(meas<0)=0;
