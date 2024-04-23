function ELEMENT = get_C_matrix( MODEL, ELEMENT, MATERIAL, KINEMATICS, formul_check )

         % % --------------- FUNCTION INFO ---------------- % %

% get_C_matrix builds the constituve law matrix associated to the chosen
% model. To do so Voigt notation is exploited. 
%
%        ELEMENT = get_C_matrix
%                   ( MODEL, ELEMENT, MATERIAL, KINEMATICS, formul_check )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% ELEMENT             [struct]      ELEMENT structure               [multi]
% MATERIAL            [struct]      MATERIAL structure              [multi]
% KINEMATICS          [struct]      KINEMATICS structure            [multi]
% formul_check        [char]        'UL' or 'TL'                    [-]
%
% -------------------------------------------------------------------------
% Output arguments:
% ELEMENT               [struct]    struct containing parameters for     
%                                   the ELEMENT of the structure    [multi]
%
% -------------------------------------------------------------------------

% --- Variable initialization
E = MATERIAL.E;
nu = MATERIAL.nu;
mu = E/(2*(1+nu));
J = KINEMATICS.J;
F = KINEMATICS.Def_gradient;
[ ~ ,dim] = size(MODEL.nodes);

if strcmpi(MATERIAL.type,'svk') == 1 

    if strcmpi(formul_check,'ul') == 1 

        ELEMENT.c = (E/(1-nu^2)).*[ 1  nu 0; 
                                nu 1  0;
                                0  0  (1-nu)/2];

        % push forward for c matrix
        C = zeros(dim,dim,dim,dim);

        C(1,1,1,1) = ELEMENT.c(1,1);
        C(1,1,2,2) = ELEMENT.c(1,2);
%         C(1,1,1,2) = ELEMENT.c(1,3);
        C(2,2,1,1) = ELEMENT.c(2,1);
        C(2,2,2,2) = ELEMENT.c(2,2);
%         C(2,2,1,2) = ELEMENT.c(2,3);
%         C(1,2,1,1) = ELEMENT.c(3,1);
%         C(1,2,2,2) = ELEMENT.c(3,2);
        C(1,2,1,2) = ELEMENT.c(3,3);
%         C(1,1,2,1) = C(1,1,1,2);
%         C(2,2,2,1) = C(2,2,1,2);
%         C(2,1,1,1) = C(1,2,1,1);
%         C(2,1,2,2) = C(1,2,2,2);
        C(2,1,1,2) = C(1,2,1,2);
        C(2,1,2,1) = C(1,2,1,2);
        C(1,2,2,1) = C(1,2,1,2);
        
        for l = 1:2
            for m = 1:2
                for n = 1:2
                    for p = 1:2
                        c(l,m,n,p) = 1/J * (F(l,1) * F(m,1) * F(n,1) * F(p,1) * C(1,1,1,1) ...
                            + F(l,1) * F(m,1) * F(n,2) * F(p,2) * C(1,1,2,2) ...
                            + F(l,2) * F(m,2) * F(n,1) * F(p,1) * C(2,2,1,1) ...
                            + F(l,2) * F(m,2) * F(n,2) * F(p,2) * C(2,2,2,2) ...
                            + F(l,1) * F(m,2) * F(n,1) * F(p,2) * C(1,2,1,2) ...
                            + F(l,2) * F(m,1) * F(n,1) * F(p,2) * C(2,1,1,2) ...
                            + F(l,2) * F(m,1) * F(n,2) * F(p,1) * C(2,1,2,1) ...
                            + F(l,1) * F(m,2) * F(n,2) * F(p,1) * C(1,2,2,1)); ...
        %                     + F(l,1) * F(m,1) * F(n,1) * F(p,2) * C(1,1,1,2) ...
        %                     + F(l,2) * F(m,2) * F(n,1) * F(p,2) * C(2,2,1,2) ...
        %                     + F(l,1) * F(m,2) * F(n,1) * F(p,1) * C(1,2,1,1) ...
        %                     + F(l,1) * F(m,2) * F(n,2) * F(p,2) * C(1,2,2,2) ...
        %                     + F(l,1) * F(m,1) * F(n,2) * F(p,1) * C(1,1,2,1) ...
        %                     + F(l,2) * F(m,2) * F(n,2) * F(p,1) * C(2,2,2,1) ...
        %                     + F(l,2) * F(m,1) * F(n,1) * F(p,1) * C(2,1,1,1) ...
        %                     + F(l,2) * F(m,1) * F(n,2) * F(p,2) * C(2,1,2,2));
                    end
                end
            end
        end
                            
        ELEMENT.c_voigt = [c(1,1,1,1), c(1,1,2,2), c(1,1,1,2)
                           c(2,2,1,1), c(2,2,2,2), c(2,2,1,2)
                           c(1,2,1,1), c(1,2,2,2), c(1,2,1,2)];

    elseif strcmpi(formul_check,'tl') == 1 
        
        ELEMENT.c = (E/(1-nu^2)).*[ 1  nu 0; 
                                nu 1  0;
                                0  0  (1-nu)/2];
    
        ELEMENT.c_voigt = ELEMENT.c;
    
    else

        fprintf('\n%s is not a valid solution type, specify one between:\n - UL (Updated Lagrangian) \n - TL (Total Lagrangian)\n', formul_check)
        error('Solution type is not valid.')
    end

elseif strcmpi(MATERIAL.type,'neohookean') == 1 

    if strcmpi(formul_check,'ul') == 1 

        lambda_prime=2*mu/(J^2);
        mu_prime = mu/(J^2);
        ELEMENT.t = MATERIAL.thickness/J;
        ELEMENT.c = lambda_prime*kron(eye(2),eye(2))+2*mu_prime*eye(4);
        ELEMENT.c_voigt = [lambda_prime+2*mu_prime, lambda_prime,           0;
                           lambda_prime,            lambda_prime+2*mu_prime,0;
                           0,                       0,                      mu_prime];

    elseif strcmpi(formul_check,'tl') == 1 

        fprintf('\nThe neohookean material model is not available with a Total Lagrangian formulation, please change the value to "UL"\n')
        error('Formulation impossible to use')

    else
       
        fprintf('\n%s is not a valid solution type, specify one between:\n - UL (Updated Lagrangian) \n - TL (Total Lagrangian)\n', formul_check)
        error('Solution type is not valid.')

    end

else
    
    fprintf('\n%s is not a valid material type, specify one between:\n - SVK (Saint Venant-Kirchhoff) \n - NeoHookean \n', MATERIAL.type)
    error('Material type type is not valid.')
end

