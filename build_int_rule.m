function [int_rule] = build_int_rule(n_integration_pts)

int_rule = struct();

switch(n_integration_pts)
    case 1
        int_rule.x = [0];
        int_rule.w = [2];
    case 2
        int_rule.x = [-1/sqrt(3) 1/sqrt(3)];
        int_rule.w = [1 1];
    case 3
        int_rule.x = [-sqrt(0.6) 0 sqrt(0.6)];
        int_rule.w = [5/9 8/9 5/9];
    otherwise
        disp('This number of integration points cannot be used');
end
