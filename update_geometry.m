function MODEL = update_geometry( MODEL )

% --- Reshape U_unc
dU_res = reshape(MODEL.dU,2,length(MODEL.dU)/2)';

% --- Update geometry for step i-th
MODEL.xy = MODEL.xy+dU_res;

% --- Save total displacements
MODEL.U = MODEL.U + MODEL.dU;

% --- Reorder displacements
[q,r] = size(MODEL.nodes);    % q = number of nodes, r = coordinates of the nodes
MODEL.UxUy = reshape(MODEL.U, r, q)';

