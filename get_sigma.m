function ELEMENT = get_sigma( ELEMENT, MATERIAL, KINEMATICS, formul_check, el_count, ip_count )
% Return: Cauchy stresses (needed here to evaluate the  geometric 
% stiffness)

% Note: evaluation of stresses and elasticity tensor depend on the
% material model 

% cambiare sigma con "stress" ????
% fare push forward ELEMENT.c
% svk vale per small strains

if strcmpi(MATERIAL.type,'svk') == 1 || strcmpi(MATERIAL.type,'plane_stress_linear') == 1

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
        % trovo strain green lagrange
        ELEMENT.E = 0.5.* (F'*F - eye(size(F)));
        % Voigt notation
        ELEMENT.E_voigt = [ELEMENT.E(1,1); ELEMENT.E(2,2); 2*ELEMENT.E(1,2)];
        % get 2nd PK stress tensor
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

        % same formulation, tbh works only for ul
        mu = MATERIAL.E/(2*(1+MATERIAL.nu));
        b = KINEMATICS.Def_gradient*KINEMATICS.Def_gradient';
        ELEMENT.sigma_tens = mu * (b - eye(size(ELEMENT.sigma_tens))/(KINEMATICS.J^2));
        ELEMENT.sigma_voigt = [ELEMENT.sigma_tens(1,1); ELEMENT.sigma_tens(2,2); ELEMENT.sigma_tens(1,2)];
        ELEMENT.sigmas(el_count,ip_count).sigma = ELEMENT.sigma_tens;

    else
       
        error('%s is not a valid solution type, specify one between:\n - UL (Updated Lagrangian) \n - TL (Total Lagrangian)\n', formul_check)
    
    end

else
    
    error('%s is not a valid material type, specify one between:\n - SVK (for Saint Venant-Kirchhoff) \n - NeoHookean \n', MATERIAL.type)

end

