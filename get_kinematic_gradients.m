function [KINEMATICS, JAC_det] = get_kinematic_gradients(KINEMATICS, ELEMENT, xi_input, eta_input, formul_check, nl_check, el_count, ip_count )

% --- Computational shape functions derivatives evaluation
KINEMATICS.dN_eval = ELEMENT.shape_function_derivatives(xi_input,eta_input);

% --- Jacobians 
KINEMATICS.Jacobian = ELEMENT.XY_nodes' * KINEMATICS.dN_eval;  % [ MATERIAL ]
KINEMATICS.jacobian = ELEMENT.xy_nodes' * KINEMATICS.dN_eval;  % [ spatial  ]

% --- Jacobian determinant
KINEMATICS.detJ = det(KINEMATICS.Jacobian); % [ MATERIAL ]
KINEMATICS.detj = det(KINEMATICS.jacobian); % [ spatial  ]

% --- Inverse jacobian
Jacobian_inverse = KINEMATICS.Jacobian\eye(size(KINEMATICS.Jacobian));
jacobian_inverse = KINEMATICS.jacobian\eye(size(KINEMATICS.jacobian));

% --- Physical shape functions derivatives
KINEMATICS.dN_material = (KINEMATICS.dN_eval * Jacobian_inverse);   % [ MATERIAL ]
KINEMATICS.dN_spatial = (KINEMATICS.dN_eval * jacobian_inverse);    % [ spatial  ]

% --- Deformation gradient
KINEMATICS.Def_gradient = ELEMENT.xy_nodes' * KINEMATICS.dN_material;

% --- Deformation gradient determinant ( J )
KINEMATICS.J = det(KINEMATICS.Def_gradient);

if strcmpi(nl_check,'nonlinear') == 1   
    % Save Deformation gradient for all element in all integration points
    KINEMATICS.F(el_count,ip_count).F_mat = KINEMATICS.Def_gradient;
end

% Assign only the right determinant as output
if strcmpi(formul_check,'ul') == 1 
    JAC_det = KINEMATICS.detj;
else
    JAC_det = KINEMATICS.detJ;
end

end