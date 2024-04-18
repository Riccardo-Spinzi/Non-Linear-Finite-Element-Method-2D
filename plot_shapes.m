function plot_shapes( MODEL )

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

    plot_x_coord_un = [plot_x_coord; plot_x_coord(1)]; % altrimenti non disegna l'ultima linea dell'elemento, potei cercare qualche function migliore di plot
    plot_y_coord_un = [plot_y_coord; plot_y_coord(1)]; % idem per y

    % undeformed
    subplot(1,2,1)
    plot(plot_x_coord_un, plot_y_coord_un, '.-b')

    plot_x_coord_def = plot_x_coord + U_def_x;
    plot_y_coord_def = plot_y_coord + U_def_y;

    plot_x_coord_def = [plot_x_coord_def; plot_x_coord_def(1)]; % altrimenti non disegna l'ultima linea dell'elemento, potei cercare qualche function migliore di plot
    plot_y_coord_def = [plot_y_coord_def; plot_y_coord_def(1)]; % idem per y

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







