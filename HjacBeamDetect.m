function H = HjacBeamDetect(mu)
%   HjacBeamDetect(mu) - outputs the linearized measurement jacobian
%   INPUTS
%       mu      [x y teta x1 y1 x2 y2 ... xi yi]'
%       
% 
%   OUTPUTS
%       H     [ddxi/dx ddxi/dy ddxi/dteta ddxi/dx1 ddxi/dy1 ...; 
%              ddyi/dx ddyi/dy ddyi/dteta ddyi/dx1 ddyi/dy1 ...; ]
% 
%   Cornell University
%   Autonomous Mobile Robots
%   Homework #6
%   Scher, Guy

Nm = (length(mu) - 3)/2;
theta = mu(3);
H = zeros(Nm*2,3+Nm*2);
for i=1:2:Nm*2
    H(i,3+i)   =  cos(theta); %ddx/dxi
    H(i,4+i)   =  sin(theta); %ddx/dyi
    H(i+1,3+i) = -sin(theta); %ddy/dxi
    H(i+1,4+i) =  cos(theta); %ddy/dyi   
    H(i, 3)    = -(mu(3+i)-mu(1))*sin(theta)+(mu(4+i)-mu(2))*cos(theta); %ddx/dteta
    H(i+1, 3)  = -(mu(3+i)-mu(1))*cos(theta)-(mu(4+i)-mu(2))*sin(theta); %ddy/dteta
end
% delta_x
H(1:2:end,1) = -cos(theta); %ddx/dxp
H(1:2:end,2) = -sin(theta); %ddx/dyp
% delta_y
H(2:2:end,1) =  sin(theta);
H(2:2:end,2) = -cos(theta);
