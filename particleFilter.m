function [X1] = particleFilter(X0, u1, Z1, R, Q0, g_full, h_full)
% adapted particleFilter: output the new set of particles given previous set 
%                         for the FAST-SLAM 1.0 problem
%
%   INPUTS
%       X0           previous set of particles  (& weights)
%       u1           current command
%       Z1           current measurement vector
%       R            state model noise
%       Q0           measurement noise (will take care of NaN in function)
%       g_full       function pointer to system dynamics inputs: (x,u)
%       h_full       function pointer to measurement model. inputs: (x)
%
%   OUTPUTS
%       X0           new set of particles (& weights)
%
%   Cornell University
%   Autonomous Mobile Robots
%   Homework 6
%   Scher, Guy

x = X0.x; w = ones(size(X0.w)); Sig = X0.S;
[Ns, Np] = size(x);
Nm = (Ns - 3)/2; % # of features
x_bar = zeros(size(X0.x));
Sig_bar = zeros(size(Sig));
X1 = X0;

Q = Q0; % measurement noise, take only valid (not nan)

Hjac   = @(x) HjacBeamDetect(x); % pointer to the linearized measurement model

for i=1:Np
    % prediction
    x_tmp = g_full(x(1:3,i), u1) + sqrt(R(1:3,1:3))*randn(3,1);  %integrateOdom and sample for pose
    x_bar(1:3,i) = x_tmp(:,end);
    for j=1:2:(Nm*2)
        if(isnan(Z1(j)) || isnan(Z1(j+1)) )
            % not seen, don't update
            x_bar(3+j,i)     = x(3+j,i);  %xj
            x_bar(3+j+1,i)   = x(3+j+1,i); % yj
            Sig_bar(2*j-1:2*j+2,i) = Sig(2*j-1:2*j+2,i); % sigma
        else
            % we have measurement. First time?
            if( all( Sig(2*j-1:2*j+2,i) == 0)  )
                xy = x_bar(1:2,i); teta = x_bar(3,i);
                T = [cos(teta) -sin(teta); sin(teta) cos(teta)];
                xy_j = T*Z1(j:j+1) + xy; %inverse measurement model
                x_bar(3+j,i)     = xy_j(1);  %xj initialize feature
                x_bar(3+j+1,i)   = xy_j(2); % yj
                H = Hjac([x_bar(1:3,i); x_bar(3+j:3+j+1,i)]); %H(x,mj)
                H = H(:, 4:5); % how they relate to xj and yj
                sig_j = (H \ Q(j:j+1,j:j+1)) * inv(H)'; % inv(H)
                Sig_bar(2*j-1:2*j+2,i) = reshape(sig_j, 4,1);
                w(i) = w(i)*1; % we don't add yet because we still don't know how valid it is
            else
                % not first time, update
                prev_m = x(3+j:3+j+1,i);
                prev_sig = reshape(Sig(2*j-1:2*j+2,i),2,2);
                H = Hjac([x_bar(1:3,i); prev_m]); %H(x,mj)
                H = H(:, 4:5);
                QQ = H*prev_sig*H' + Q(j:j+1,j:j+1);
                K = prev_sig*H' / QQ;
                exp_meas = h_full([x_bar(1:3,i); prev_m]); %expected meas
                xy_j = prev_m + K*(Z1(j:j+1)-exp_meas); %
                x_bar(3+j,i)     = xy_j(1);  %xj initialize feature
                x_bar(3+j+1,i)   = xy_j(2); % yj
                sig_j = (eye(2)-K*H)*prev_sig;
                Sig_bar(2*j-1:2*j+2,i) = reshape(sig_j, 4,1);
                tmp = QQ\(exp_meas-Z1(j:j+1));
                w(i) = w(i) * 1/(2*pi*det(QQ)) * exp(-(exp_meas-Z1(j:j+1))'*tmp/2 ); %Multivariate normal probability density function, like mvnpdf
            end
        end
    end
end
% Normalize to form a probability distribution (sum = 1).
w_sum = sum(w);
if(w_sum == 0)
    X1.w = ones(1, Np)/Np; %we're in a pickle, so re-intialize to make it even again
else
    X1.w = w./sum(w);
end
% resampling:
v = rand(Np,1); 
wc = cumsum(X1.w'); 
[~ ,ind1] = sort([v; wc]);
ind=find(ind1<=Np)-(0:Np-1)'; 
X1.x(:,:) = x_bar(:,ind);
X1.S(:,:) = Sig_bar(:,ind);

