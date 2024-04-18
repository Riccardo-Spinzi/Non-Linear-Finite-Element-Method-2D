function [ MODEL, SOL, POST ] = linear_solver( MODEL, SOL, POST, MATERIAL, ELEMENT )
 
% ---Load vector (eval external forces)
[ MODEL, SOL ] = build_fext( MODEL, SOL );
 
% ---Stiffness matrix
MODEL = build_K_matrix( MODEL, MATERIAL, ELEMENT, SOL.type, 'ul' );
 
% ---Solve linear sistem
MODEL = solve_structure( MODEL, SOL.type );

% ---Get stresses 
MODEL = stress_recovery ( MODEL, SOL, POST, MATERIAL, ELEMENT );