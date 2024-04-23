function Display_results( MODEL, SOL, POST )

         % % --------------- FUNCTION INFO ---------------- % %

% Display_results outputs the useful results to the command window, so that
% scrolling to search them is no longer a burden.
%
%                   Display_results( MODEL, SOL, POST )
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% SOL                 [struct]      SOL structure                   [multi]
% POST                [struct]      POST structure                  [multi]
% -------------------------------------------------------------------------
% Output arguments:
%
% -------------------------------------------------------------------------

[m,n] = size(MODEL.elements);
fprintf('---- RESULTS ----\n')
fprintf('Problem type: %s\n', upper(SOL.type))
fprintf('Number of elements: %d\n', m)
fprintf('Number of nodes per element: %d\n', n)
fprintf('-----------------\n')
if strcmpi(SOL.type,'linear')
    fprintf('Maximum absolute Displacement in x direction: %6.4f\n', max(abs(MODEL.UxUy(:,1))))
    fprintf('Maximum absolute Displacement in y direction: %6.4f\n', max(abs(MODEL.UxUy(:,2))))
    fprintf('Maximum Absolute sigma_xx: %6.4f\n', max(max(abs(MODEL.sigma_xx))))
    fprintf('Maximum Absolute sigma_yy: %6.4f\n', max(max(abs(MODEL.sigma_yy))))
    fprintf('Maximum Absolute sigma_xy: %6.4f\n', max(max(abs(MODEL.sigma_xy))))
    fprintf('Strain energy U calculated: %6.4f\n', MODEL.Strain_energy)
else
    fprintf('Results referred to last load step:\n')
    fprintf('Maximum absolute Displacement in x direction: %6.4f\n', max(abs(POST.STEP(end).Ux)))
    fprintf('Maximum absolute Displacement in y direction: %6.4f\n', max(abs(POST.STEP(end).Uy)))
    fprintf('Maximum Absolute sigma_xx: %6.4f\n', max(max(abs(MODEL.sigma_xx))))
    fprintf('Maximum Absolute sigma_yy: %6.4f\n', max(max(abs(MODEL.sigma_yy))))
    fprintf('Maximum Absolute sigma_xy: %6.4f\n', max(max(abs(MODEL.sigma_xy))))
    fprintf('Strain energy U calculated: %6.4f\n', MODEL.Strain_energy)
end
fprintf('---------------- NOTE:\n')
fprintf(['Units of measure depend on the input which was given initially:\n' ...
    'for example, if the position was given in mm and the forces on the\n' ...
    'nodes in N, the displacements will be in mm, the stresses in MPa \n' ...
    'and the strain energy in Nmm, always check for consistency in units\n' ...
    'of measure\n'])



