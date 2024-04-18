function MODEL = set_pointers( MODEL )

% ogni riga rappresenta un elemento e ogni colonna i dof associati a
% quell'elemento, uso il vettore pos per aiutarmi.
[m,n] = size(MODEL.pos);                
[q,r] = size(MODEL.elements);
MODEL.ptrs = zeros(q,MODEL.eltype*n);

% order nodes in a vector following the order declaration of the elements
x = reshape(MODEL.elements',1,q*r);  
% get a vector with as columns the dofs of the different nodes
y = MODEL.pos(x,:)';
% reorder these dofs as a row vector
z = reshape(y,1,q*r*n);
% use them to fill the ptrs matrix
MODEL.ptrs = reshape(z,r*n,q)';
