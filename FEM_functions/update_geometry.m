function MODEL = update_geometry( MODEL )

         % % --------------- FUNCTION INFO ---------------- % %

% update_geometry uses the solution obtained at the current loading step
% (only non lineare case) to update the position of the nodes (spatial
% configuration) and the total displacement of the structure.
% 
%                   MODEL = update_geometry( MODEL )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
%
% -------------------------------------------------------------------------
% Output arguments:
% MODEL               [struct]      struct containing parameters     
%                                   for the MODEL of the structure  [multi]
% -------------------------------------------------------------------------

% --- Reshape U_unc
dU_res = reshape(MODEL.dU,2,length(MODEL.dU)/2)';

% --- Update geometry for step i-th
MODEL.xy = MODEL.xy+dU_res;

% --- Save total displacements
MODEL.U = MODEL.U + MODEL.dU;

% --- Reorder displacements
[q,r] = size(MODEL.nodes);    % q = number of nodes, r = coordinates of the nodes
MODEL.UxUy = reshape(MODEL.U, r, q)';

