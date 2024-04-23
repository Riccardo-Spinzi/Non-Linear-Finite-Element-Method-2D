function [ MODEL, SOL ] = build_fext( MODEL, SOL )

         % % --------------- FUNCTION INFO ---------------- % %

% build_fext takes as input the force shape vector Fhat and gives as output
% the force vector to apply to the structure. It is the whole force vector
% if the problem is linear, otherwise it multiplies the shape force vector
% with the current lambda parameter to deflate the total load and create a
% load step, applying only partially the force shape factor MODEL.Fhat.
%
%             [ MODEL, SOL ] = build_fext( MODEL, SOL )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% SOL                 [struct]      SOL structure                   [multi]
%
% -------------------------------------------------------------------------
% Output arguments:
% MODEL               [struct]      struct containing parameters     
%                                   for the MODEL of the structure  [multi]
% SOL                 [struct]      struct containing parameters     
%                                   for the SOL of the structure    [multi]
% -------------------------------------------------------------------------
 
switch lower(SOL.type)
    case 'linear'
        MODEL.F = MODEL.Fhat;
    case 'nonlinear'
        SOL.lambda = SOL.lambda + SOL.dlambda;
        MODEL.F = MODEL.Fhat*(SOL.lambda);
    otherwise
        error('Unknown method.');
end