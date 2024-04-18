function MODEL = stress_recovery( MODEL, SOL, POST, MATERIAL, ELEMENT )

Barlow = MODEL.Barlow_points;
nl_check = SOL.type;
formul_check = SOL.formul;

switch Barlow 
    case 0
        nint_p = length( MODEL.int_rule.x );
        int_rule = MODEL.int_rule;
    case 1
        nint_p = length( MODEL.int_rule.x ) - Barlow;
        [int_rule] = build_int_rule(nint_p);
    otherwise
        error(['The specified value for Barlow points integration is not valid, specify:\n- 1 to activate Barlow points integration' ...
            '\n- 0 to avoid Barlow points integration'])
end

sigma_xx_plot = zeros(MODEL.nels,nint_p^2);

[m,n] = size(MODEL.elements);
X = zeros(m,n);
Y  = zeros(m,n);

for i = 1 : MODEL.nels

    count = 1;

    % ---Matrix of nodal coordinates (spatial & material)
    el_nodes = MODEL.elements( i, : );

    ELEMENT.XY_nodes= [ MODEL.XY( el_nodes, 1 ) MODEL.XY( el_nodes, 2 ) ];
    ELEMENT.xy_nodes= [ MODEL.xy( el_nodes, 1 ) MODEL.xy( el_nodes, 2 ) ];

    for iG = 1 : nint_p

        for jG = 1 : nint_p

            xi = int_rule.x(jG);
            eta = int_rule.x(iG);
            nodes = MODEL.elements(i,:);

            % --- Get kinematic gradients
            [KINEMATICS, ~] = get_kinematic_gradients( [], ELEMENT, xi, eta, formul_check, nl_check, i, count );       
       
            % --- Get elasticity tensor
            ELEMENT = get_C_matrix( ELEMENT, MATERIAL, KINEMATICS, formul_check );        
    
            % --- Get B matrix
            [ Bc, ~ ] = get_B_matrix( KINEMATICS, nl_check ); 

            % --- Recover useful displacements
            U = MODEL.UxUy(nodes,:);
            [r,c] = size(U);
            U_calc = reshape(U',r*c,1);

            % --- Recover stress at integration points
            stress_voigt = ELEMENT.c * Bc * U_calc; 
            sigma_vect_xx(count) = stress_voigt(1);
            sigma_vect_yy(count) = stress_voigt(2);
            sigma_vect_xy(count) = stress_voigt(3);

            count = count + 1;        

        end

    end

    sigma_xx_plot(i,:) = sigma_vect_xx;
    sigma_yy_plot(i,:) = sigma_vect_yy;
    sigma_xy_plot(i,:) = sigma_vect_xy;


end

if n == 8 && Barlow == 0
    
    sigma_xx_plot(:,[2,4,5,6,8]) = [];
    sigma_yy_plot(:,[2,4,5,6,8]) = [];
    sigma_xy_plot(:,[2,4,5,6,8]) = [];

end

% Manipulate to get the right distribution for PATCH
% sigma_xx_plot(:,[2,3]) = sigma_xx_plot(:,[3,2]);
sigma_xx_plot(:,[3,4]) = sigma_xx_plot(:,[4,3]);
% sigma_yy_plot(:,[2,3]) = sigma_yy_plot(:,[3,2]);
sigma_yy_plot(:,[3,4]) = sigma_yy_plot(:,[4,3]);
% sigma_xy_plot(:,[2,3]) = sigma_xy_plot(:,[3,2]);
sigma_xy_plot(:,[3,4]) = sigma_xy_plot(:,[4,3]);

MODEL.sigma_xx = sigma_xx_plot;
MODEL.sigma_yy = sigma_yy_plot;
MODEL.sigma_xy = sigma_xy_plot;
