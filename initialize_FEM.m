function [ MODEL, SOL ] = initialize_FEM( MODEL, SOL )

% ---Initial coordinates
MODEL.XY = MODEL.nodes;
MODEL.xy = MODEL.XY;

if strcmp(SOL.type,'nonlinear') == 1
    % ---Load parameters
    SOL.lambda = 0;
    SOL.lambdas= [];
end