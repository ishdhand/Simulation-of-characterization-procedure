%% Unitary Constraint
% Constraints the matrix to be unitary by pre and post
% multiplying with appropriate real diagnol matrices
%% Inputs:
%       relativeAmplitudesMatrix: mxm matrix which records the amplitude of the
%           transition matrix relative to the elements of the first row and
%           column. The entries of the first row and column are 1.
%       phasesMatrix: mxm matrix with the phases of the
%           transition matrix. The first row and column are 0.
%% outputs:
%       TransitionMatrix: mxm unitary matrix
%% Procedure

function TransitionMatrix = UnitaryConstraint(relativeAmplitudesMatrix,phasesMatrix)

% size of interferometer
m = length(relativeAmplitudesMatrix);

% the matrix M_\mu
Mu = relativeAmplitudesMatrix.*exp(1i*phasesMatrix);

% Solve for the sigmas
B1 = [1; zeros(m-1,1)];
sigmaVector = linsolve(Mu,B1);

% matrix S
S = diag(sqrt(abs(sigmaVector)));

% solve for the rhos
B2 = B1/sigmaVector(1);
rhoVector = linsolve(Mu',B2);

% set the first element manually to 1.
rhoVector(1) = 1;

% matrix R
R = diag(sqrt(abs(rhoVector)));

% the form of the transition matrix
TransitionMatrix = R*Mu*S;

end

