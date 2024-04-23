function plot_shapes( MODEL )

         % % --------------- FUNCTION INFO ---------------- % %

% plot_shapes displays the displacements undergoing in the structure at the
% final step of integration. Since it shows qualitative results, the 
% contributions of the nodes which are not at the corners of an element are
% killed.
% 
%                       plot_shapes( MODEL )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
%
% -------------------------------------------------------------------------
% Output arguments:
%
% -------------------------------------------------------------------------

% --- Declare useful variables
elements = MODEL.elements;
nodes = MODEL.nodes;
UxUy = MODEL.UxUy;
n_visual = 0.75;
[m,n] = size(MODEL.elements); % m = number of elements, n = number of nodes

% --- Plot displacements 
figure(1)
subplot(1,2,1)
hold on
for i = 1 : m

    nodes_el = MODEL.elements(i,:);
    plot_x_coord = nodes(nodes_el,1);
    plot_y_coord = nodes(nodes_el,2);
    
    U_def_x = n_visual*UxUy(nodes_el,1);
    U_def_y = n_visual*UxUy(nodes_el,2);

    if length(nodes) > 4 % just take corners to plot the figure
        plot_x_coord(5:end) = [];
        plot_y_coord(5:end) = [];
        U_def_x(5:end) = [];
        U_def_y(5:end) = [];
    end

    plot_x_coord_un = [plot_x_coord; plot_x_coord(1)]; 
    plot_y_coord_un = [plot_y_coord; plot_y_coord(1)]; 

    % undeformed
    subplot(1,2,1)
    plot(plot_x_coord_un, plot_y_coord_un, '.-b')

    plot_x_coord_def = plot_x_coord + U_def_x;
    plot_y_coord_def = plot_y_coord + U_def_y;

    plot_x_coord_def = [plot_x_coord_def; plot_x_coord_def(1)]; 
    plot_y_coord_def = [plot_y_coord_def; plot_y_coord_def(1)];

    % deformed
    subplot(1,2,1)
    plot(plot_x_coord_def, plot_y_coord_def, 'r')
    
end

% --- Patch plot
if n > 4
    elements(:,5:end) = [];
end
figure(1)
subplot(1,2,2)
col = sqrt(UxUy(:,1).^2 + UxUy(:,2).^2);
patch('Faces', elements,'Vertices',nodes + UxUy, ...
    'FaceVertexCData', col, 'FaceColor', 'interp')
colorbar

sgtitle('Displacement plots') 






