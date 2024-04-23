function [ MODEL, SOL ] = initialize_FEM( MODEL, SOL )

         % % --------------- FUNCTION INFO ---------------- % %

% initialize_FEM sets the parameters at the start of the simulation, such
% are the position of the nodes at the initial step and the lambda factor
% for the non linear solution procedure.
%
%             [ MODEL, SOL ] = initialize_FEM( MODEL, SOL )
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

% ---Initial coordinates
MODEL.XY = MODEL.nodes;
MODEL.xy = MODEL.XY;

if strcmp(SOL.type,'nonlinear') == 1
    % ---Load parameters
    SOL.lambda = 0;
    SOL.lambdas= [];
end