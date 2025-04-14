function [ MODEL, POST, MATERIAL, SOL, ELEMENT ] = set_model( INPUT, NR_method )

         % % --------------- FUNCTION INFO ---------------- % %

% set_model is the only function that is not thought to solve the problem,
% but rather to set the model in a way that it can be interpreted by the
% program and such that the problem can be solved in the most organised way
% possibile.
%
%  [ MODEL, POST, MATERIAL, SOL, ELEMENT ] = set_model( INPUT, NR_method )
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
% POST                [struct]      struct containing parameters     
%                                   for the POST of the structure   [multi]
% MATERIAL            [struct]      struct containing parameters for     
%                                   the MATERIAL of the structure   [multi]
% SOL                 [struct]      struct containing parameters     
%                                   for the SOL of the structure    [multi]
% ELEMENT             [struct]      struct containing parameters for   
%                                   the ELEMENT of the structure    [multi]
% -------------------------------------------------------------------------

% struct initialization
MODEL = struct();
POST = struct();
MATERIAL = struct();
SOL = struct();
ELEMENT = struct();

% --- MODEL 
[m,n] = size(INPUT.elements);       % m is the number of elements, n is the number of nodes per element
[q,~] = size(INPUT.load);           % q is the number of loads applied, ~ will be always 3 for program syntax
[s,l] = size(INPUT.nodes);          % s is the number of nodes, l is the number of dofs of the single node
[r,~] = size(INPUT.spc);            % r is the number of constraints
MODEL.ndof = s*l;
MODEL.nfree_dofs = MODEL.ndof - length(INPUT.spc);
MODEL.nels = m;
MODEL.eltype = n;
MODEL.elements = INPUT.elements;
MODEL.nodes = INPUT.nodes;
MODEL.pos = (reshape(1:MODEL.ndof,l,length(MODEL.nodes)))';     
for i = 1 : r
    MODEL.constr_dofs(i,1) = MODEL.pos(INPUT.spc(i,1),INPUT.spc(i,2));
end
MODEL.free_dofs = (1:MODEL.ndof)';                              
MODEL.free_dofs(MODEL.constr_dofs) = [];
MODEL.int_rule = build_int_rule(INPUT.integration_pts);         
MODEL.K = zeros(MODEL.ndof,MODEL.ndof);
MODEL.F = zeros(MODEL.ndof,1);                                  
MODEL.Fhat = zeros(MODEL.ndof,1);
for i = 1 : q
    MODEL.Fhat(MODEL.pos(INPUT.load(i,1),INPUT.load(i,2))) = INPUT.load(i,3);
end
MODEL.K_cons = zeros(MODEL.nfree_dofs,MODEL.nfree_dofs);     
MODEL.F_cons = zeros(MODEL.nfree_dofs,1);                    
MODEL.U = zeros(MODEL.ndof,1);                                 

if strcmp(INPUT.sol_type,'nonlinear') == 1
    MODEL.fint = zeros(MODEL.ndof,1);
    MODEL.res = zeros(MODEL.ndof,1);
    MODEL.res_cons = zeros(MODEL.nfree_dofs,1);                
    MODEL.dU = zeros(MODEL.ndof,1);
    MODEL.dU_cons = zeros(MODEL.nfree_dofs,1);                 
end

MODEL.NR_method = NR_method;
MODEL.Barlow_points = INPUT.Barlow_points;

% --- MATERIAL
MATERIAL.type = INPUT.mat_type;
MATERIAL.E = INPUT.mat_prop(1);
MATERIAL.nu = INPUT.mat_prop(2);
MATERIAL.thickness = INPUT.mat_prop(3);

% --- SOL 
SOL.type = INPUT.sol_type;
SOL.lambda_max = INPUT.NR_param(1);
SOL.dlambda = INPUT.NR_param(2);
SOL.nincr = [];
SOL.formul = INPUT.formul;                             
if strcmp(INPUT.sol_type,'nonlinear')==1
    SOL.norm_res_max = INPUT.NR_param(3);
    SOL.nitermax = INPUT.NR_param(4);
    SOL.lambda = 0;
    SOL.lambdas = SOL.lambda:SOL.dlambda:SOL.lambda_max;
end

% --- POST 
if strcmp(INPUT.sol_type,'nonlinear')==1
    for i = 1 : length(SOL.lambdas) - 1               
        POST.STEP(i).X = zeros(MODEL.nels,length(MODEL.nodes));
        POST.STEP(i).Y = zeros(MODEL.nels,length(MODEL.nodes));
        POST.STEP(i).x = zeros(MODEL.nels,length(MODEL.nodes));       
        POST.STEP(i).y = zeros(MODEL.nels,length(MODEL.nodes));
        POST.STEP(i).Ux = zeros(MODEL.nels,length(MODEL.nodes));      
        POST.STEP(i).Uy = zeros(MODEL.nels,length(MODEL.nodes));
    end
end


% --- ELEMENT
switch(n)
    case 4 
        ELEMENT.shape_functions = @(xi,eta) [0.25.*(1-xi).*(1-eta);
                                             0.25.*(1+xi).*(1-eta);
                                             0.25.*(1+xi).*(1+eta);
                                             0.25.*(1-xi).*(1+eta);];

        ELEMENT.shape_function_derivatives = @(xi,eta) [-0.25.*(1-eta), -0.25.*(1-xi);
                                                        0.25.*(1-eta), -0.25.*(1+xi);
                                                        0.25.*(1+eta), 0.25.*(1+xi);
                                                        -0.25.*(1+eta), 0.25.*(1-xi)];
    case 8
        ELEMENT.shape_functions = @(xi,eta) [0.25.*(1-xi).*(1-eta).*(-xi-eta-1);
                                             0.25.*(1+xi).*(1-eta).*(xi-eta-1);
                                             0.25.*(1+xi).*(1+eta).*(xi+eta-1);
                                             0.25.*(1-xi).*(1+eta).*(-xi+eta-1);
                                             0.5.*(1-xi.^2).*(1-eta);
                                             0.5.*(1+xi).*(1-eta.^2);
                                             0.5.*(1-xi.^2).*(1+eta);
                                             0.5.*(1-xi).*(1-eta.^2);];

        ELEMENT.shape_function_derivatives = @(xi,eta) [-((eta + 2*xi)*(eta - 1))/4, -((2*eta + xi)*(xi - 1))/4;
                                                        ((eta - 2*xi)*(eta - 1))/4, ((xi + 1)*(2*eta - xi))/4;
                                                        ((eta + 2*xi)*(eta + 1))/4, ((2*eta + xi)*(xi + 1))/4;
                                                        -((eta - 2*xi)*(eta + 1))/4, -((xi - 1)*(2*eta - xi))/4;
                                                        xi*(eta - 1), xi^2/2 - 1/2;
                                                        1/2 - eta^2/2, -eta*(xi + 1);
                                                        -xi*(eta + 1), 1/2 - xi^2/2;
                                                        eta^2/2 - 1/2,eta*(xi - 1)];

    otherwise
        disp(['il numero di nodi per ogni elemento non è corretto, inserirne...' ...
            ' 4 per elementi bilineari, 8 per elementi biquadratici ']);
end
ELEMENT.order = n;
ELEMENT.XY_nodes = zeros(ELEMENT.order,l);          
ELEMENT.xy_nodes = zeros(ELEMENT.order,l);          
ELEMENT.c = [];
ELEMENT.sigma_voigt = zeros(3,1);                   
ELEMENT.sigma_tens = zeros(2,2);
ELEMENT.t = MATERIAL.thickness;
if l ~= 2 
    error('Only 2D elements can be implemented up to now, check INPUT.nodes')
end







