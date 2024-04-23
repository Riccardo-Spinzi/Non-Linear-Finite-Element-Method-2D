function ELEMENT = get_sigma( ELEMENT, MATERIAL, KINEMATICS, formul_check, el_count, ip_count )

         % % --------------- FUNCTION INFO ---------------- % %

% get_sigma evaluates the stresses at the integration points, to later use
% to build the tangent stiffness matrix. 
%
%            ELEMENT = get_sigma
%                         ( ELEMENT, MATERIAL, KINEMATICS, 
%                           formul_check, el_count, ip_count )
%
% -------------------------------------------------------------------------
% Input arguments: 
% ELEMENT             [struct]      ELEMENT structure               [multi]
% MATERIAL            [struct]      MATERIAL structure              [multi]
% KINEMATICS          [struct]      KINEMATICS structure            [multi]
% formul_check        [char]        'UL' or 'TL'                    [-]
% el_count            [1x1 double]  element number                  [-]
% ip_count            [1x1 double]  int. point number               [-]
%
% -------------------------------------------------------------------------
% Output arguments:
% ELEMENT             [struct]      struct containing parameters for   
%                                   the ELEMENT of the structure    [multi]
%
% -------------------------------------------------------------------------

if strcmpi(MATERIAL.type,'svk') == 1 

    if strcmpi(formul_check,'ul') == 1  

        F = KINEMATICS.Def_gradient;
        F_inv = F\eye(size(F));
        ELEMENT.e = 0.5*(eye(size(F))-F_inv'*F_inv);
        ELEMENT.e_voigt = [ELEMENT.e(1,1); ELEMENT.e(2,2); 2*ELEMENT.e(1,2)];
        ELEMENT.sigma_voigt = ELEMENT.c_voigt*ELEMENT.e_voigt;  
        ELEMENT.sigma_tens = [ELEMENT.sigma_voigt(1),   ELEMENT.sigma_voigt(3);
                              ELEMENT.sigma_voigt(3),   ELEMENT.sigma_voigt(2)];
        ELEMENT.sigmas(el_count,ip_count).sigma = ELEMENT.sigma_tens;
        
    elseif strcmpi(formul_check,'tl') == 1 
        
        F = KINEMATICS.Def_gradient;
        ELEMENT.E = 0.5.* (F'*F - eye(size(F)));
        ELEMENT.E_voigt = [ELEMENT.E(1,1); ELEMENT.E(2,2); 2*ELEMENT.E(1,2)];
        ELEMENT.sigma_voigt = ELEMENT.c_voigt*ELEMENT.E_voigt;  
        ELEMENT.sigma_tens = [ELEMENT.sigma_voigt(1),   ELEMENT.sigma_voigt(3);
                              ELEMENT.sigma_voigt(3),   ELEMENT.sigma_voigt(2)];
        ELEMENT.sigmas(el_count,ip_count).sigma = ELEMENT.sigma_tens;

    else

        error('%s is not a valid solution type, specify one between:\n - UL (Updated Lagrangian) \n - TL (Total Lagrangian)\n', formul_check)
    
    end

elseif strcmpi(MATERIAL.type,'neohookean') == 1 

    if strcmpi(formul_check,'ul') == 1 

        mu = MATERIAL.E/(2*(1+MATERIAL.nu));
        b = KINEMATICS.Def_gradient*KINEMATICS.Def_gradient';
        ELEMENT.sigma_tens = mu * (b - eye(size(ELEMENT.sigma_tens))/(KINEMATICS.J^2));
        ELEMENT.sigma_voigt = [ELEMENT.sigma_tens(1,1); ELEMENT.sigma_tens(2,2); ELEMENT.sigma_tens(1,2)];
        ELEMENT.sigmas(el_count,ip_count).sigma = ELEMENT.sigma_tens;

    elseif strcmpi(formul_check,'tl') == 1 

        fprintf('\nThe neohookean material model is not available with a Total Lagrangian formulation, please change the value to "UL"\n')
        error('Formulation impossible to use')

    else
       
        fprintf('%s is not a valid solution type, specify one between:\n - UL (Updated Lagrangian) \n - TL (Total Lagrangian)\n', formul_check)
        error('Formulation not found')

    end

else
    
    frpintf('%s is not a valid material type, specify one between:\n - SVK (for Saint Venant-Kirchhoff) \n - NeoHookean \n', MATERIAL.type)
    error('Material model not found')

end

