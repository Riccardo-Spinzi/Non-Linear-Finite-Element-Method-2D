function MODEL = build_K_matrix( MODEL, MATERIAL, ELEMENT, nl_check, formul_check )

% ---Integration rule parameters
xGauss      = MODEL.int_rule.x;
wGauss      = MODEL.int_rule.w;
nint_p      = length( xGauss );
[ ~, dim ]  = size(MODEL.nodes); % dim is useful to generalize K's dimensions ( for 3D elements )

for i = 1 : MODEL.nels

    % ---Initialize
    Kc_el= zeros( dim * MODEL.eltype, dim * MODEL.eltype );     

    % ---Matrix of nodal coordinates (spatial & material)
    el_nodes = MODEL.elements( i, : );

    ELEMENT.XY_nodes= [ MODEL.XY( el_nodes, 1 ) MODEL.XY( el_nodes, 2 ) ];
    ELEMENT.xy_nodes= [ MODEL.xy( el_nodes, 1 ) MODEL.xy( el_nodes, 2 ) ];

    % ---Integration point number ( 1 to 4 ) 
    count = 1;

    % ---Integrate stiffness matrix
    for iG= 1 : nint_p                  % cicle on x integration points
    
        for jG= 1 : nint_p              % cicle on y integration points
    
            xi  = xGauss( jG);          
    
            eta = xGauss( iG);          
    
            % --- Get kinematic gradients
            [KINEMATICS, ~] = get_kinematic_gradients( [], ELEMENT, xi, eta, formul_check, nl_check, i, count );       
           
            % --- Get elasticity tensor
            ELEMENT = get_C_matrix( ELEMENT, MATERIAL, KINEMATICS, formul_check );        
    
            % --- Get B matrix
            [ Bc, ~ ] = get_B_matrix( KINEMATICS, nl_check ); 

            % --- Save B matrix for later use 
%             MODEL.B_mat(i).Bc(iG,jG).B = Bc;

            % --- Coefficient for integration in the computational domain
            dOm = wGauss( iG) * wGauss( jG ) * KINEMATICS.detJ * ELEMENT.t;
    
            % --- Stiffness matrix 
            Kc_el = Kc_el + ( Bc' * ELEMENT.c* Bc) * dOm;
    
            count = count + 1;
        end
    
    end

    % ---Assemble contribution
    ptrs= MODEL.ptrs( i, : );
    MODEL.K( ptrs, ptrs) = MODEL.K( ptrs, ptrs ) + Kc_el;
    
   
end

% --- Save c matrix
MODEL.c = ELEMENT.c;

% ... Save data only for free dofs (constrained structure)
constr_dofs = MODEL.constr_dofs;
MODEL.K_cons = MODEL.K;
MODEL.F_cons = MODEL.F;                     

% ---Impose constraints
MODEL.K_cons( constr_dofs, : ) = [];
MODEL.K_cons( :, constr_dofs) = [];
MODEL.F_cons( constr_dofs, : ) = [];            

