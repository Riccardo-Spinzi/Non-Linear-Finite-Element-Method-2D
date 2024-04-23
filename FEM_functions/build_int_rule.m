function [int_rule] = build_int_rule(n_integration_pts)

         % % --------------- FUNCTION INFO ---------------- % %

% build_int_rule gets as input the number of integration points and gives
% as output the parameters that will be used to integrate the problem, such
% as position and weights of the integration points.
%
%             [int_rule] = build_int_rule(n_integration_pts)
%
% -------------------------------------------------------------------------
% Input arguments: 
% n_integration_pts   [1 x 1 double] Number of integration points   [-]
%
% -------------------------------------------------------------------------
% Output arguments:
% int_rule            [struct]       Structure with position and 
%                                    weights of integration points  [multi]
% -------------------------------------------------------------------------
 
int_rule = struct();

switch(n_integration_pts)
    case 1
        int_rule.x = 0;
        int_rule.w = 2;
    case 2
        int_rule.x = [-1/sqrt(3) 1/sqrt(3)];
        int_rule.w = [1 1];
    case 3
        int_rule.x = [-sqrt(0.6) 0 sqrt(0.6)];
        int_rule.w = [5/9 8/9 5/9];
    otherwise
        fprintf('\nThis number of integration points cannot be used\n Please choose a number between 1 and 3.\n');
        error('Impossible to continue the simulation')
end
