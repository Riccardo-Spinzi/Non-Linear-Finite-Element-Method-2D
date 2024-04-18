function [ MODEL, POST, MATERIAL, SOL, ELEMENT ] = set_model( INPUT, NR_method )

% struct initialization
MODEL = struct();
POST = struct();
MATERIAL = struct();
SOL = struct();
ELEMENT = struct();

% --- MODEL 
[m,n] = size(INPUT.elements);       % m è il numero di elementi, n il numero di nodi associati ad un elemento
[q,~] = size(INPUT.load);           % q è il numero di forze, r sarà sempre = 3 per sintassi dell'algoritmo
[s,l] = size(INPUT.nodes);          % s è il numero di nodi, l è il numero di gdl del nodo
MODEL.ndof = s*l;
MODEL.nfree_dofs = MODEL.ndof - length(INPUT.spc);
MODEL.nels = m;
MODEL.eltype = n;
MODEL.elements = INPUT.elements;
MODEL.nodes = INPUT.nodes;
MODEL.pos = (reshape(1:MODEL.ndof,l,length(MODEL.nodes)))';     
for i = 1 : length(INPUT.spc)
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
    MODEL.res_cons = zeros(MODEL.nfree_dofs,1);                % qui potrei mettere gia il conto dei gdl <--- che vuol dire LOL
    MODEL.dU = zeros(MODEL.ndof,1);
    MODEL.dU_cons = zeros(MODEL.nfree_dofs,1);                 % qui potrei mettere gia il conto dei gdl
end

MODEL.NR_method = NR_method;
MODEL.Barlow_points = INPUT.Barlow_points;

% --- MATERIAL
MATERIAL.type = INPUT.mat_type;
MATERIAL.E = INPUT.E;
MATERIAL.nu = INPUT.nu;
MATERIAL.thickness = INPUT.t;

% --- SOL 
SOL.type = INPUT.sol_type;
SOL.lambda_max = INPUT.lambda_max;
SOL.dlambda = INPUT.dlambda;
SOL.nincr = [];
SOL.formul = INPUT.formul;                              % choose between UL and TL
if strcmp(INPUT.sol_type,'nonlinear')==1
    SOL.norm_res_max = INPUT.norm_res_max;
    SOL.nitermax = INPUT.nitermax;
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







