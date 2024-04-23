function [ MODEL, SOL ] = linear_solver( MODEL, SOL, MATERIAL, ELEMENT )

         % % --------------- FUNCTION INFO ---------------- % %

% linear_solver computes the solution of the linear problem, evaluating 
% the stiffness matrix and inverting it to solve the problem.
%
%    [ MODEL, SOL ] = linear_solver( MODEL, SOL, MATERIAL, ELEMENT )
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
% MODEL               [struct]      struct containing parameters     
%                                   for the MODEL of the structure  [multi]
% SOL                 [struct]      struct containing parameters     
%                                   for the SOL of the structure    [multi]
% -------------------------------------------------------------------------

% ---Load vector (eval external forces)
[ MODEL, SOL ] = build_fext( MODEL, SOL );
 
% ---Stiffness matrix
MODEL = build_K_matrix( MODEL, MATERIAL, ELEMENT, SOL.type, 'ul' );
 
% ---Solve linear sistem
MODEL = solve_structure( MODEL, SOL.type );

% ---Get stresses 
MODEL = stress_recovery ( MODEL, SOL, MATERIAL, ELEMENT );

% Get strain energy
MODEL.Strain_energy = 0.5 * MODEL.U' * MODEL.K * MODEL.U;




