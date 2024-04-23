function MODEL = solve_structure( MODEL, nl_check )

         % % --------------- FUNCTION INFO ---------------- % %

% solve_structure inverts the stiffness matrix and solves the linear
% problem to get displacements.  
%
%                MODEL = solve_structure( MODEL, nl_check )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% nl_check            [char]        'linear' or 'nonlinear'         [-]
% -------------------------------------------------------------------------
% Output arguments:
% MODEL               [struct]      struct containing parameters     
%                                   for the MODEL of the structure  [multi]
% -------------------------------------------------------------------------

if strcmpi(nl_check,'linear') == 1 
    % --- Solve the problem
    MODEL.U_cons = MODEL.K_cons \ MODEL.F_cons;
    
    % --- Expand displacements to the global vector
    MODEL.U = zeros( MODEL.ndof, 1);
    MODEL.U( MODEL.free_dofs ) = MODEL.U_cons;
    
    % --- Reorder displacements
    [q,r] = size(MODEL.nodes);    % q = number of nodes, r = coordinates of the nodes
    MODEL.UxUy = reshape(MODEL.U, r, q)'; 
    
else
    % --- Solve the problem
    MODEL.dU_cons = MODEL.K_cons \ MODEL.res_cons;
    
    % --- Expand displacements to the global vector
    MODEL.dU = zeros( MODEL.ndof, 1);
    MODEL.dU( MODEL.free_dofs ) = MODEL.dU_cons;
    
end
