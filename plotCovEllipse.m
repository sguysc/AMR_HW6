function objHands = plotCovEllipse(Mu,P,Nsig,plotOpts,figHandle)
% Plots the N-sigma ellipse of a 2-dimensional Gaussian distribution
% REQUIRED INPUTS:
%    Mu       - [2x1] Mean vector (or [1x2])
%    P        - [2x2] Covariance matrix
%
% OPTIONAL INPUTS:
%    Nsig     - Vector of confidence intervals to plot (default = 1)
%               Example: 
%               Nsig = [1,2,3] plots the 1,2,and 3-sigma ellipses
%               corresponding to 68.2%, 95.4%, and 99.7% confidence 
%               intervals, respectively
%
%   plotOpts  - [1x2*N] Cell array of 'N' plot options 
%               (default = matlabs default) 
%               Example:  
%               plotOpts = [{'color'},{'b'},{'linestyle'},{'-.'},{'linewidth'},{2}];
%               plots each ellipse with the following options:
%               color = blue, linestyle = dash-dot, linewidth = 2  
%
%   figHandle - Figure handle to plot into (default = current figure)
%             
% OUTPUTS:
%   objHands - handles of objects plotted
%
% Written by: Kevin Wyffels (klw94@cornell.edu), 4/1/2013

numPts = 100;                % Number of points to discretize ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check input arguments and assign defaults where necessary 
if ~exist('figHandle','var') || isempty(figHandle), figHandle = gcf; end
if ~exist('Nsig','var') || isempty(Nsig), Nsig = 1; end
    
[r,c] = size(P);
dim = length(Mu);
if (r~=2)||(c~=2), error('Covariance matrix must be [2x2]'), end
if (dim ~= 2), error('Mean vector must be of length 2'), end
  
%% Compute Ellipse Parameters
% discretize 360 degrees into 'numPts' points
theta = 0:2*pi/(numPts-1):2*pi;  

% Compute the ellipse parameters from the covariance matrix
[e1 e2 phi] = computeEllipseParams(P);

%% Plot Ellipses
% Pre-allocate memory for plot object handles
objHands = zeros(1,length(Nsig));

% Scale major and minor axes by Nsig
a = Nsig*e1;
b = Nsig*e2;

figure(figHandle); hold on,
for n = 1:length(Nsig)
    % Evaluate ellipse at discretized angles
    x = Mu(1) + a(n)*cos(theta)*cos(phi) - b(n)*sin(theta)*sin(phi);
    y = Mu(2) + a(n)*cos(theta)*sin(phi) + b(n)*sin(theta)*cos(phi);

    % Plot ellipse
    objHands(n) = plot(x,y);
end
hold off

if exist('plotOpts','var')
    % Set plot options
    Nopts = numel(plotOpts)/2;
    plotOpts = reshape(plotOpts,2,Nopts)';
    for n = 1:Nopts,
        set(objHands,cell2mat(plotOpts(n,1)),cell2mat(plotOpts(n,2)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                - Define Member Functions -                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e1 e2 phi] = computeEllipseParams(P)
% Computes the parameters of a 1-sigma ellipse centered around the origin
% given a covariance matrix
% Inputs:
%        P - 2x2 covariance matrix
% Outputs:
%       e1    - Major radius of covariance ellipse
%       e2    - Minor radius of covariance ellipse
%       theta - [radians] counter-clockwise angle from positive x-axis 
%               to ellipse major axis (e1) 

% Compute eigenvalues and eigenvectors of covariance matrix
[evec,eval] = eig(P);

% Distinguish between major and minor axes.
[eval,ix] = sort(diag(eval),'descend');
evec = evec(:,ix);

% Compute angle between positive x-axis and major axis
e1 = sqrt(eval(1));
e2 = sqrt(eval(2));
phi = atan2(evec(2,1),evec(1,1));