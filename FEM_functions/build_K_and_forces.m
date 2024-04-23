function [MODEL, ELEMENT, KINEMATICS] = build_K_and_forces( MODEL, MATERIAL, ELEMENT, SOL, niter )

         % % --------------- FUNCTION INFO ---------------- % %

% build_K_and_forces is the core of the non linear solver. It takes as
% input a large variety of parameters and assembles the tanget stiffness 
% matrix of the structure, in order to solve the non linear problem at 
% every NR iteration. Also, it evaluates stresses through the constitutive
% law specified by the user in the INPUT struct. After the tangent
% stiffness matrix has been built, the program imposes the constraints and
% decides which stiffness matrix to use depending on the NR method chosen
% at the start (in the FEM_RUNME.m file).
%
%  [MODEL, ELEMENT, KINEMATICS] = build_K_and_forces
%                                  ( MODEL, MATERIAL, ELEMENT, SOL, niter )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% MATERIAL            [struct]      MATERIAL structure              [multi]
% ELEMENT             [struct]      ELEMENT structure               [multi]
% SOL                 [struct]      SOL structure                   [multi]
% iter                [1x1 double]  Number of NR iterations         [-]
% -------------------------------------------------------------------------
% Output arguments:
% MODEL               [struct]      struct containing parameters for    
%                                   the MODEL of the structure      [multi]
% ELEMENT             [struct]      struct containing parameters for     
%                                   the ELEMENT of the structure    [multi]
% KINEMATICS          [struct]      struct containing parameters for     
%                                   the KINEMATICS of the structure [multi]
% -------------------------------------------------------------------------

% --- Init
MODEL.K = zeros( MODEL.ndof, MODEL.ndof );
MODEL.fint = zeros( MODEL.ndof, 1 );

% --- Integration rule parameters
xGauss = MODEL.int_rule.x;
wGauss = MODEL.int_rule.w;
nint_p = length( xGauss );
KINEMATICS = struct();

% --- Other parameters
nl_check = SOL.type;
formul_check = SOL.formul;
nincr = SOL.nincr;

for i = 1 : MODEL.nels

    % --- Initialize
    Kc_el = zeros( 2 * MODEL.eltype, 2 * MODEL.eltype );
    Kg_el = zeros( 2 * MODEL.eltype, 2 * MODEL.eltype );
    fint = zeros( 2 * MODEL.eltype, 1 );

    % --- Matrix of nodal coordinates (spatial & material)
    el_nodes = MODEL.elements( i, : );
    ELEMENT.XY_nodes = [ MODEL.XY( el_nodes, 1 ) MODEL.XY( el_nodes, 2 ) ];
    ELEMENT.xy_nodes = [ MODEL.xy( el_nodes, 1 ) MODEL.xy( el_nodes, 2 ) ];

    count = 1;

    % --- Integrate stiffness matrix
        for iG = 1 : nint_p

            for jG = 1 : nint_p

                xi = xGauss( jG );
                eta = xGauss( iG );

                % --- Get kinematic gradients
                [KINEMATICS, JAC_det] = get_kinematic_gradients(KINEMATICS, ELEMENT, xi, eta, formul_check, nl_check, i, count);

                % --- Get elasticity tensor
                ELEMENT = get_C_matrix( MODEL, ELEMENT, MATERIAL, KINEMATICS, formul_check );  
                ELEMENT = get_sigma( ELEMENT, MATERIAL, KINEMATICS, formul_check, i, count);
                
                % --- Get B matrix
                [ Bc, Bg ] = get_B_matrix( KINEMATICS, nl_check, formul_check );

                % --- Coefficient for integration in the comput domain
                dOm = wGauss( iG ) * wGauss( jG ) * JAC_det * ELEMENT.t;

                % --- Stiffness matrix
                Kc_el = Kc_el + ( Bc' * ELEMENT.c_voigt * Bc ) * dOm; 
                Kg_el = Kg_el + kron( Bg' * ELEMENT.sigma_tens * Bg, eye( size(ELEMENT.sigma_tens) ) ) * dOm;

                % --- Internal forces
                fint = fint + Bc' * ELEMENT.sigma_voigt * dOm; 
                
                % aggiorno count
                count = count + 1;
            end

        end

    % --- Assemble contribution
    ptrs = MODEL.ptrs( i, : );
    MODEL.K( ptrs, ptrs ) = MODEL.K( ptrs, ptrs ) + Kc_el + Kg_el;
    MODEL.fint( ptrs, 1 ) = MODEL.fint( ptrs, 1 ) + fint;

end

% ... Residual
MODEL.res = MODEL.F - MODEL.fint;

% ... Save data only for free dofs (constrained structure)
constr_dofs = MODEL.constr_dofs;
MODEL.F_cons = MODEL.F;
MODEL.res_cons = MODEL.res;

switch lower(MODEL.NR_method)
    case 'nr'
         MODEL.K_cons = MODEL.K;
    case 'nr_mod'
        if niter == 1
            MODEL.K_NR_mod = MODEL.K;
        end
        MODEL.K_cons = MODEL.K_NR_mod;
    case 'nr_initial'
        if nincr == 1 && niter == 1
            MODEL.K_NR_initial = MODEL.K;
        end
        MODEL.K_cons = MODEL.K_NR_initial;
    otherwise
        fprintf(['\n Invalid tangent stiffness update method. Choose between:\n' ...
            '- NR for Newton-Rhapson method (update K at every iteration)\n' ...
            '- NR_mod for modified Newton-Rhapson method (update K at the beginning of every load step)\n' ...
            '- NR_initial for initial tangent stiffness method (Evaluate K only once at the start)\n'])
        error('NR method specified not found.')
end

% --- Impose constraints
MODEL.K_cons( constr_dofs, : ) = [];
MODEL.K_cons( :, constr_dofs ) = [];
MODEL.res_cons( constr_dofs, : ) = [];
MODEL.F_cons( constr_dofs, : ) = [];
