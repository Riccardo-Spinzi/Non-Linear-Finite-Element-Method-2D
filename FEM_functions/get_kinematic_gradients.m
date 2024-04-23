function [KINEMATICS, JAC_det] = get_kinematic_gradients(KINEMATICS, ELEMENT, xi_input, eta_input, formul_check, nl_check, el_count, ip_count )

         % % --------------- FUNCTION INFO ---------------- % %

% get_kinematic_grandients is another core part of the program, where
% derivatives of the shape functions, material and spatial jacobians and 
% deformation gradient are evaluated. To do so, a lot of parameters are 
% needed as input.
%
%  [KINEMATICS, JAC_det] = get_kinematic_gradients
%                           (KINEMATICS, ELEMENT, xi_input, eta_input, 
%                            formul_check, nl_check, el_count, ip_count )
%
% -------------------------------------------------------------------------
% Input arguments: 
% ELEMENT             [struct]      ELEMENT structure               [multi]
% KINEMATICS          [struct]      KINEMATICS structure            [multi]
% xi_input            [1x1 double]  x coordinate of the int. point  [-]
% eta_input           [1x1 double]  y coordinate of the int. point  [-]
% formul_check        [char]        'UL' or 'TL'                    [-]
% nl_check            [char]        'linear' or 'nonlinear'         [-]
% el_count            [1x1 double]  element number                  [-]
% ip_count            [1x1 double]  int. point number               [-]
%
% -------------------------------------------------------------------------
% Output arguments:
% KINEMATICS          [struct]      struct containing parameters for     
%                                   the KINEMATICS of the structure [multi]
% JAC_det             [1x1 double]  determinant of the jacobian     [-]
%
% -------------------------------------------------------------------------

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