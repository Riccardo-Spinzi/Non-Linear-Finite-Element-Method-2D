function MODEL = set_pointers( MODEL )

         % % --------------- FUNCTION INFO ---------------- % %

% set_pointers uses the pos vector to associate the right pointers at every
% element, so to build correctly the tangent stiffness matrix later.
%
%                    MODEL = set_pointers( MODEL )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
%
% -------------------------------------------------------------------------
% Output arguments:
% MODEL               [struct]      struct containing parameters     
%                                   for the MODEL of the structure  [multi]
%
% -------------------------------------------------------------------------

[~,n] = size(MODEL.pos);        % n is the number of dofs associated to a node    
[q,r] = size(MODEL.elements);   % q is the number of elements, r is the number of nodes associated to one element
MODEL.ptrs = zeros(q,MODEL.eltype*n);

% order nodes in a vector following the order declaration of the elements
x = reshape(MODEL.elements',1,q*r);  
% get a vector with as columns the dofs of the different nodes
y = MODEL.pos(x,:)';
% reorder these dofs as a row vector
z = reshape(y,1,q*r*n);
% use them to fill the ptrs matrix
MODEL.ptrs = reshape(z,r*n,q)';
