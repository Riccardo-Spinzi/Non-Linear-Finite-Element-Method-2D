function [ MODEL, SOL, POST ] = fem_solver( INPUT, NR_method )

         % % --------------- FUNCTION INFO ---------------- % %

% fem_solver manages the workflow of the program, being the higher function
% in the hierarchy. From there, data is stored, controlled and sent where 
% it needs to be in order to solve the structural problem. 
%
%          [ MODEL, SOL, POST ] = fem_solver( INPUT, NR_method )
%
% -------------------------------------------------------------------------
% Input arguments: 
% INPUT               [struct]      INPUT structure                 [multi]
% NR_method           [char]        kind of NR method used to assign
%                                   the tangent stiffness matrix    [-]
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

% ---Set model
 
[ MODEL, POST, MATERIAL, SOL, ELEMENT ] = set_model( INPUT, NR_method );
 
% ---Set pointers
 
MODEL = set_pointers( MODEL );
 
 % ---Initialize
 
[ MODEL, SOL ] = initialize_FEM( MODEL, SOL );
 
% ---Solution type
 
switch SOL.type
   case 'linear'
           [ MODEL, SOL ] = linear_solver( MODEL, SOL, MATERIAL, ELEMENT );
   case 'nonlinear'
           [ MODEL, SOL, POST ] = nonlinear_solver( MODEL, SOL, POST, MATERIAL, ELEMENT );
   otherwise
           fprintf('\nUnknown method.\n Please choose between "linear" or "nonlinear"\n');
           error('Unknow problem type.')
end
