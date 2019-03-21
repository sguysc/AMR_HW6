function RunSLAM(rec_file)
% RunSLAM(rec_file): outputs two plots for the EKF-SLAM and the FAST-SLAM 1.0
%
%   INPUTS
%       rec_file     name of input file in format:
%                    time, V, W, measx1, measy2, ...measxn, measyn
%
%   OUTPUTS
%       plots
%
%   Cornell University
%   Autonomous Mobile Robots
%   Homework 6
%   Scher, Guy

    % colors for plots
    C = [   0         0    1.0000; ...
        1.0000         0         0; ...
            0    1.0000         0; ...
            0         0    0.1724; ...
        1.0000    0.1034    0.7241; ...
        1.0000    0.8276         0; ...
            0    0.3448         0; ...
        0.5172    0.5172    1.0000; ...
        0.6207    0.3103    0.2759; ...
            0    1.0000    0.7586];
    if nargin == 0
       rec_file = 'hw6SLAM.txt'; 
    end
    data = importdata(rec_file);
    
    Nm = (size(data,2) - 3)/2; %-3 pose, and each one has x,y 
    
    t = data(:,1);
    controls = data(:, 2:3);
    feats = data(:, 4:end);

    N_particles = 50;

    %% EFK SLAM
    mu = [zeros(3,1); zeros(Nm*2, 1)]'; %x,y,theta,m1x,m1y,m2x ...
    deadReck = zeros(3,1)';
    %covariance of the full matrix
    Sig = 1e2*eye(3+Nm*2); Sig(1,1)=0.01; Sig(2,2)=0.01; Sig(3,3)=0.01;
    Sig = reshape(Sig, 1, (3+Nm*2)^2);
    
    R = [0.01 0.001 0.002; 0.001 0.01 0.002; 0.002 0.002 0.015];
    %combine the process noise of the full matrix
    R = [R zeros(3, Nm*2); zeros(Nm*2, 3) 0.001*eye(Nm*2)];
    Q = 0.1*eye(Nm*2); % x,y ~N(0,0.1) measurements

    for i=1:length(t)
        dt = 0;
        if(i>1)
            dt = t(i)-t(i-1);
        end
        u1 = dt*controls(i, :)'; % convert V,W to d,phi
        Sig0 = reshape(Sig(i, :),3+Nm*2,3+Nm*2);
        mu0 = mu(i,:)';
        newpose = integrateOdom(deadReck(i,:)', u1);
        deadReck = [deadReck; newpose(:,end)'];
        % the measurement
        Z1     = feats(i, :)';
        % dirty fix to make invalid meas. is NaN
            Z1(Z1 == 0) = NaN;
        g_full = @(x,u) integrateOdom(x, u); %pointer to full dynamics
        Gjac   = @(x,u) GjacDiffDrive(x, u); %pointer to linearized dynamics
        h_full = @(x) hBeamDetect(x); % pointer to the measurement model
        Hjac   = @(x) HjacBeamDetect(x); % pointer to the linearized measurement model
        % where the estimation magic happens
        [mu1, Sig1] = extendedKalmanFilter(mu0, u1, Sig0, Z1, R, Q, ...
                                       g_full, Gjac, h_full, Hjac);
        %store the estimation
        mu    = [mu;  mu1'];                          
        Sig = [Sig; reshape(Sig1,1,(3+Nm*2)^2)]; %store it back as 

        if i==length(t)
            figHandle=figure(1); p=zeros(9,1);
            p(1)=plot(mu(:,1), mu(:,2), ':*', 'Color', C(1,:),'DisplayName', 'Path SLAM'); hold on;
            p(2)=plot(mu(end,1), mu(end,2), 'o', 'MarkerSize', 18, 'LineWidth', 3, 'Color', C(2,:),'DisplayName', 'Cur. pose'); 
            p(3)=plot([mu(end,1) mu(end,1)+0.5*cos(mu(end,3))], ...
                 [mu(end,2) mu(end,2)+0.5*sin(mu(end,3))], '-', 'Color', C(3,:), 'MarkerSize', 18, 'LineWidth', 3,'DisplayName', 'Heading'); 
            p(4)=plot(mu(1,1), mu(1,2), 'x', 'MarkerSize', 18, 'Color', C(4,:),'DisplayName', 'Start pose'); 
            p(5)=plot(deadReck(:,1), deadReck(:,2), ':s', 'Color', C(5,:),'DisplayName', 'DeadReck');
            plotOpts = [{'color'},{C(6,:)},{'linestyle'},{':'},{'linewidth'},{2},{'DisplayName'}, {'1\sigma'}];
            p(6) = plotCovEllipse([mu(end,1), mu(end,2)], ...
                        Sig1(1:2,1:2),1,plotOpts,figHandle); hold on;
            for j=4:4:2*Nm
                p(7)=plot(mu(end,j), mu(end,j+1), '<', 'Color', C(7,:),'DisplayName', 'a marker'); 
                text(mu(end,j), mu(end,j+1), num2str((j-2)/2)); 
                p(8)=plot(mu(end,j+2), mu(end,j+3), '>', 'Color', C(8,:),'DisplayName', 'b marker'); 
                text(mu(end,j+2), mu(end,j+3), num2str((j)/2)); 
                p(9)=plot([mu(end,j) mu(end,j+2)], [mu(end,j+1) mu(end,j+3)], '-', 'LineWidth', 2, 'Color', C(9,:),'DisplayName', 'Wall'); 
                text((mu(end,j)+mu(end,j+2))/2, (mu(end,j+1)+mu(end,j+3))/2, ...
                    ['wall' num2str(j/4)]); 
                if(Sig1(j)~=Sig(1,4))
                    objHands = plotCovEllipse([mu(end,j), mu(end,j+1)], ...
                        Sig1(j:j+1,j:j+1),1,plotOpts,figHandle); hold on;
                end
                if(Sig1(j+2)~=Sig(1,4))
                    objHands = plotCovEllipse([mu(end,j+2), mu(end,j+3)], ...
                        Sig1(j+2:j+3,j+2:j+3),1,plotOpts,figHandle); hold on;
                end
            end
            grid on; hold off; axis equal; title(['EKF-SLAM: Step #' num2str(i)]);
            xlabel('X [m]'); ylabel('Y [m]'); lgd = legend(p); legend('Location','best'); title(lgd,'Legend');
    %         pause
        end
    end
    %% FAST-SLAM
    mu=[];
    
    particles.x = zeros(3+Nm*2, N_particles); %xp,yp,tetap, x1,y1, ...
    particles.S = zeros(Nm*4, N_particles); %sig11,sig12,sig21,sig22, ....
    particles.w = ones(1, N_particles)/N_particles;
    for i=1:length(t)
        dt = 0;
        if(i>1)
            dt = t(i)-t(i-1);
        end
        u1 = dt*controls(i, :)'; % convert V,W to d,phi
        % the measurement
        Z1     = feats(i, :)';
        % dirty fix to make invalid meas. is NaN
        Z1(Z1 == 0) = NaN;
        g_full = @(x,u) integrateOdom(x, u); %pointer to full dynamics
        h_full = @(x) hBeamDetect(x); % pointer to the measurement model
        X1 = particleFilter(particles(i), u1, Z1, R, Q, g_full, h_full);
        particles(i+1) = X1; % expand the list
        mu = [mu; mean(particles(i+1).x(:, :),2)'];

        if i==length(t) %for plotting last step
            figHandle=figure(2); p=zeros(8,1);
            p(1)=plot(mu(:,1), mu(:,2), ':*', 'Color', C(1,:),'DisplayName', 'Path Fast-SLAM'); hold on;
            p(2)=plot(mu(end,1), mu(end,2), 'o', 'MarkerSize', 18, 'LineWidth', 3, 'Color', C(2,:),'DisplayName', 'Cur. pose'); 
            p(3)=plot([mu(end,1) mu(end,1)+0.5*cos(mu(end,3))], ...
                 [mu(end,2) mu(end,2)+0.5*sin(mu(end,3))], '-', 'Color', C(3,:), 'MarkerSize', 18, 'LineWidth', 3,'DisplayName', 'Heading'); 
            p(4)=plot(mu(1,1), mu(1,2), 'x', 'MarkerSize', 18, 'Color', C(4,:),'DisplayName', 'Start pose'); 
            p(5)=plot(deadReck(1:i+1,1), deadReck(1:i+1,2), ':s', 'Color', C(5,:),'DisplayName', 'DeadReck');
            p(9)=plot(particles(i+1).x(1, :), particles(i+1).x(2, :),'.', 'Color', C(9,:), 'DisplayName', 'Particles (pose)');
            for j=4:4:2*Nm
                p(6)=plot(mu(end,j), mu(end,j+1), '<', 'Color', C(7,:),'DisplayName', 'a marker'); 
                text(mu(end,j), mu(end,j+1), num2str((j-2)/2)); 
                p(7)=plot(mu(end,j+2), mu(end,j+3), '>', 'Color', C(8,:),'DisplayName', 'b marker'); 
                text(mu(end,j+2), mu(end,j+3), num2str((j)/2)); 
                p(8)=plot([mu(end,j) mu(end,j+2)], [mu(end,j+1) mu(end,j+3)], '-', 'LineWidth', 2, 'Color', C(9,:),'DisplayName', 'Wall'); 
                text((mu(end,j)+mu(end,j+2))/2, (mu(end,j+1)+mu(end,j+3))/2, ...
                    ['wall' num2str(j/4)]); 
            end
            grid on; hold off; axis equal; title(['FAST-SLAM Step #' num2str(i)]);
            xlabel('X [m]'); ylabel('Y [m]'); lgd = legend(p); legend('Location','best'); title(lgd,'Legend'); hold off
    %         pause
        end
    end
end


