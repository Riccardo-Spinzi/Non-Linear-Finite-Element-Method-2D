function [ MODEL, SOL ] = build_fext( MODEL, SOL )

switch lower(SOL.type)
    case 'linear'
        MODEL.F = MODEL.Fhat;
    case 'nonlinear'
        SOL.lambda = SOL.lambda + SOL.dlambda;
        MODEL.F = MODEL.Fhat*(SOL.lambda);
    otherwise
        error('Unknown method.');
end