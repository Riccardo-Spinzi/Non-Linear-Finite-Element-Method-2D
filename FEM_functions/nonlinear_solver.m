function [ MODEL, SOL, POST ] = nonlinear_solver( MODEL, SOL, POST, MATERIAL, ELEMENT )

         % % --------------- FUNCTION INFO ---------------- % %

% nonlinear_solver computes the solution of the non linear problem,
% evaluating the stiffness matrices and inverting them to solve the problem
% at every load step and NR iteration. It also updates the solution at 
% every step. 
%
%     [ MODEL, SOL, POST ] = nonlinear_solver
%                              ( MODEL, SOL, POST, MATERIAL, ELEMENT )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% SOL                 [struct]      SOL structure                   [multi]
% POST                [struct]      POST structure                  [multi]
% MATERIAL            [struct]      MATERIAL structure              [multi]
% ELEMENT             [struct]      ELEMENT structure               [multi]
% 
% 
% -------------------------------------------------------------------------
% Output arguments:
% MODEL               [struct]      struct containing parameters     
%                                   for the MODEL of the structure  [multi]
% SOL                 [struct]      struct containing parameters     
%                                   for the SOL of the structure    [multi]
% POST                [struct]      struct containing parameters     
%                                   for the POST of the structure   [multi]
% -------------------------------------------------------------------------

% --- Begin incremental procedure: load steps
SOL.nincr = 1;

while ( abs( SOL.lambda_max - SOL.lambda ) > 1e-8 ) 

    norm_res = 2 * SOL.norm_res_max; % to enter NR iter
    fprintf('\nLoad step: %d of %d -- Iter: ', ...
    SOL.nincr, SOL.lambda_max / SOL.dlambda );

    % --- Apply load step (eval external forces)
    [ MODEL, SOL ] = build_fext( MODEL, SOL );

    % --- NR iterations
    niter = 1;

    % --- Tangent stiffness & internal forces
    [MODEL, ELEMENT, KINEMATICS] = build_K_and_forces( MODEL, MATERIAL, ELEMENT, SOL, niter );

    while ( norm_res > SOL.norm_res_max && niter <= SOL.nitermax )
        % --- Impose constraints and solve
        MODEL = solve_structure( MODEL, SOL.type );

        % --- Update configuration
        MODEL = update_geometry( MODEL );

        % --- Tangent stiffness & internal forces
        [MODEL, ELEMENT, KINEMATICS] = build_K_and_forces( MODEL, MATERIAL, ELEMENT, SOL, niter );
        
        % --- Check residual
        norm_res = norm( MODEL.res_cons ) / norm( MODEL.F_cons );

        % --- Update iteration
        niter = niter + 1;
    end
    fprintf(' %d \n', niter - 1);
    fprintf('Error: %e \n', norm_res );

    % --- Store results
    SOL.lambdas = [ SOL.lambdas; SOL.lambda ];
    POST = POST_update_deformed_shape( POST, MODEL, SOL.nincr );
    POST = POST_stress_recovery( MODEL, SOL, POST, MATERIAL, ELEMENT );  
    POST = POST_strain_recovery( MODEL, POST, KINEMATICS, SOL.nincr );
    
    % --- Update increment
    SOL.nincr = SOL.nincr + 1;
    
end

% Get strain energy at the final load step
MODEL.Strain_energy = 0.5 * MODEL.U' * MODEL.K * MODEL.U;





