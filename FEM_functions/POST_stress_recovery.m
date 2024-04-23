function POST = POST_stress_recovery( MODEL, SOL, POST, MATERIAL, ELEMENT )

         % % --------------- FUNCTION INFO ---------------- % %

% POST_stress_recovery gets the stesses at every integration point,
% computing again the constituive law and calling the get_sigma.m function.
%
%    POST = POST_stress_recovery( MODEL, SOL, POST, MATERIAL, ELEMENT )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% SOL                 [struct]      SOL structure                   [multi]
% POST                [struct]      POST structure                  [multi]
% MATERIAL            [struct]      MATERIAL structure              [multi]
% ELEMENT             [struct]      ELEMENT structure               [multi]
%
% -------------------------------------------------------------------------
% Output arguments:
% POST                [struct]      struct containing parameters     
%                                   for the POST of the structure   [multi]
%
% -------------------------------------------------------------------------

Barlow = MODEL.Barlow_points;
nl_check = SOL.type;
formul_check = SOL.formul;
KINEMATICS = struct();

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

% --- Other parameters
nincr = SOL.nincr;
xGauss = int_rule.x;

for i = 1 : MODEL.nels

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
                [ KINEMATICS, ~ ] = get_kinematic_gradients( KINEMATICS, ELEMENT, xi, eta, formul_check, nl_check, i, count);

                % --- Get elasticity tensor  
                ELEMENT = get_sigma( ELEMENT, MATERIAL, KINEMATICS, formul_check, i, count);
                
                % aggiorno count
                count = count + 1;
            end

        end

end

POST.STEP(nincr).stress= ELEMENT.sigmas;
