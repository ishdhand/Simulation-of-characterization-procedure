%% Remove input and output phases
% For unitary matrix M calculates an equivalent
% unitary matrix with input and output phases zero.
%% Inputs:
%       M: mxm unitary matrix
%% Outputs:
%       Mclean: mxm unitary matrix with real entries in first row and first
%       column.
%       D1: mxm diagnol matrix that was used to pre-multiply M.
%       D2: mxm diagnol matrix that was used to post-multiply M.
%% Procedure

function [Mclean,D1,D2] = removeInputOutputPhases(M)
m = length(M);
iphases = angle(M);
% phases of complex entries of M

% the phases from the first row and col will be removed
B = -[iphases(1,1:m) iphases(2:m,1)']';

% Consruct and solve a linear system to find the values of phases that need
% to be pre and post multiplied with the matrix.
% phaseeq = ...
%     [1 0 0 1 0 0;
%      1 0 0 0 1 0;
%      1 0 0 0 0 1;
%      0 1 0 1 0 0;
%      0 0 1 1 0 0;];
phaseeq = zeros(2*m-1,2*m);
phaseeq(1:m,1) = 1;
phaseeq(m+1:2*m-1,m+1) = 1;
for i=1:m
    phaseeq(i,m+i) = 1;
end
for i=1:m-1
    phaseeq(m+i,i+1) = 1;
end


% find the phases
phaseVec = linsolve(phaseeq,B);

% the input and output phase matrices
D1 = diag(exp(1i*phaseVec(1:m)));
D2 = diag(exp(1i*phaseVec(m+1:2*m)));

% matrix with input output phases removed
Mclean = D1*M*D2;

if angle(Mclean(2,2)) < 0
    Ma = abs(Mclean);
    Mp = angle(Mclean);
    Mclean = Ma.*exp(-1i*Mp);
end

end