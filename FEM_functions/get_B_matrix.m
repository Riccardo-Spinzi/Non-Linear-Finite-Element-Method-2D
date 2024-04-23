function [ Bc, Bg ] = get_B_matrix( KINEMATICS, nl_check, formul_check )

         % % --------------- FUNCTION INFO ---------------- % %

% get_B_matrix computes the matrix containing geometric and material info,
% in order to build the tangent stiffness matrix for an element. To do so,
% the shape functions derivatives are evaluated at the integration points.
% 
%     [ Bc, Bg ] = get_B_matrix( KINEMATICS, nl_check, formul_check )
%
% -------------------------------------------------------------------------
% Input arguments: 
% KINEMATICS          [struct]      KINEMATICS structure            [multi]
% nl_check            [char]        'linear' or 'nonlinear'         [-]
% formul_check        [char]        'UL' or 'TL'                    [-]
%
% -------------------------------------------------------------------------
% Output arguments:
% Bc     [3x(2*eltype) double]    Matrix with material contribution 
%                                 to the tangent stiffness matrix   [-]
% Bg     [(eltype)x2 double]      Matrix with geometric contribution 
%                                 to the tangent stiffness matrix   [-]
% -------------------------------------------------------------------------

switch lower(nl_check)
    case 'linear'
        Bc = [];
        Bg = [];

        dN=KINEMATICS.dN_spatial;

        for i = 1 : length(dN)
            Bc_temp = [dN(i,1)   0;
                  0         dN(i,2);
                  dN(i,2)   dN(i,1)];
            Bc = [Bc Bc_temp];
        end

    case 'nonlinear'

        Bc = [];
        Bg = [];
        

        % --- Check for Updated or total Lagrangian formulation
        if strcmpi(formul_check,'ul')
            
            dN=KINEMATICS.dN_spatial;
            for i = 1 : length(dN)
                Bc_temp = [dN(i,1)   0;
                      0         dN(i,2);
                      dN(i,2)   dN(i,1)];
                Bc = [Bc Bc_temp];
            end

        elseif strcmpi(formul_check,'tl')

            F = KINEMATICS.Def_gradient;
            dN=KINEMATICS.dN_material;
            
            for i = 1 : length(dN)
                Bc_temp = [F(1,1)*dN(i,1)                   F(2,1)*dN(i,1);
                           F(1,2)*dN(i,2)                   F(2,2)*dN(i,2);
                           F(1,2)*dN(i,1)+ F(1,1)*dN(i,2)   F(2,1)*dN(i,2)+ F(2,2)*dN(i,1)];
                Bc = [Bc Bc_temp];
            end

        end
        
        % Bg does not change between UL and TL
        Bg_temp = dN;
        Bg = [Bg Bg_temp']; 

    otherwise

        error('Unknown method')

end