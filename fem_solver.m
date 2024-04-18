function [ MODEL, SOL, POST ] = fem_solver( INPUT, NR_method )
    
% ---Set model
 
[ MODEL, POST, MATERIAL, SOL, ELEMENT ] = set_model( INPUT, NR_method );
 
% ---Set pointers
 
MODEL = set_pointers( MODEL );
 
 % ---Initialize
 
[ MODEL, SOL ] = initialize_FEM( MODEL, SOL );
 
% ---Solution type
 
switch SOL.type
   case 'linear'
           [ MODEL, SOL, POST ] = linear_solver( MODEL, SOL, [], MATERIAL, ELEMENT );
   case 'nonlinear'
           [ MODEL, SOL, POST ] = nonlinear_solver( MODEL, SOL, POST, MATERIAL, ELEMENT );
   otherwise
           disp('Unknown method.');
end
