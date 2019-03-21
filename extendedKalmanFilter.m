function [mu1, Sig1] = extendedKalmanFilter(mu0, u1, Sig0, Z1, R, Q0, ...
                                            g_full, Gjac, h_full, Hjac)
% extendedKalmanFilter: output the belief given previous belief
%
%   INPUTS
%       mu0          previous vector of pose
%       u1           current command
%       Sig0         previous covariance matrix
%       Z1           current measurement vector
%       R            state model noise
%       Q0           measurement noise (will take care of NaN in function)
%       g_full       function pointer to system dynamics inputs: (x,u)
%       Gjac         linearized system dynamics inputs: (x)
%       h_full       function pointer to measurement model. inputs: (x)
%       Hjac         function pointer to linearized measurements. inputs: (x)
%
%   OUTPUTS
%       mu1          current estimate of vector of pose
%       Sig1         current covariance matrix
%
%   Cornell University
%   Autonomous Mobile Robots
%   Homework 4
%   Scher, Guy

% prediction
mu_bar = g_full(mu0, u1);           %integrateOdom
mu_bar = mu_bar(:,end);             %take last measurement

G = Gjac(mu0, u1);                  %GjacDiffDrive
Sig_bar = G*Sig0*G' + R;            %advance covariance based on linearized dynamics

valid = ~isnan(Z1); % if we get NaNs then ignore that measurement
% no good measurement, continue to next one
if(isempty(valid))
    % uhhmmmm, exit with only prediction made
    mu1 = mu_bar;
    Sig1 = Sig_bar;
    disp('no valid measurements in EKF');
else
    H = Hjac(mu_bar); % set the meas. model to suite what we got
    H = H(valid,:);  %filter out the NaNs
	exp_meas = h_full(mu_bar);
    % % just test to filter out noisy meas.
%     v_tst = abs(exp_meas-Z1)>1; %exclude based on absolute error
%     % or
%     v_tst = abs(exp_meas-Z1) > (mean(exp_meas-Z1) + 3*std(exp_meas-Z1)); %exclude based on relative error
%     for i=1:length(v_tst)
%         if(v_tst(i)==1)
%             Q0(i, i) = 16; % something big as oppused to the original Q
%         end
%     end
    
    Q = Q0(valid, valid); % measurement noise, take only valid (not nan)
    % Kalman gain
    Kt = Sig_bar * H' / ( H * Sig_bar * H' + Q); % kalman gain
    % update
    mu1 = mu_bar + Kt*(Z1(valid) - exp_meas(valid));
    Sig1 = (eye(size(Sig0)) - Kt*H)*Sig_bar;
%   % for debugging and setting a breakpoint   
%     if( norm(mu1(1:2)-mu0(1:2)) > 0.5 )
%         disp('out of range');
%     end
    if(isnan(mu1) ) %isnan(rcond( H * Sig_bar * H' + Q) ))
        disp('mm');
    end
end 