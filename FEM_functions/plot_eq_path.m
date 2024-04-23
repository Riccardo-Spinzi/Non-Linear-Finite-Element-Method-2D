function plot_eq_path ( MODEL, POST, SOL)

         % % --------------- FUNCTION INFO ---------------- % %

% plot_eq_path displays the load steps during the overall integration
% process, to highlght possibile non linear behaviors of the structure.
%
%                   plot_eq_path ( MODEL, POST, SOL)
%
% -------------------------------------------------------------------------
% Input arguments: 
% MODEL               [struct]      MODEL structure                 [multi]
% POST                [struct]      POST structure                  [multi]
% SOL                 [struct]      SOL structure                   [multi]
%
% -------------------------------------------------------------------------
% Output arguments:
%
% -------------------------------------------------------------------------

switch lower(SOL.type)
    case 'linear'
        U_plot_x = [0; max(abs(MODEL.UxUy(:,1)))];
        U_plot_y = [0; max(abs(MODEL.UxUy(:,2)))];

        % plot
        figure
        subplot(1,2,1)
        plot(U_plot_x, [0; SOL.lambda_max],'-o')
        xlabel('Ux')
        ylabel('\lambda')
        grid on
        subplot(1,2,2)
        plot(U_plot_y, [0; SOL.lambda_max],'-o')
        xlabel('Uy')
        ylabel('\lambda')
        grid on

    case 'nonlinear'
        U_plot_x = zeros(size(SOL.lambdas)+ [1,0]);
        U_plot_y = zeros(size(SOL.lambdas)+ [1,0]);
        
        % get (max(abs(U))) for Ux e Uy
        for i = 1 : length(SOL.lambdas)
            U_plot_x(i+1) = max(abs(POST.STEP(i).Ux));
            U_plot_y(i+1) = max(abs(POST.STEP(i).Uy));
        end
        
        % plot
        figure
        subplot(1,2,1)
        plot(U_plot_x, [0; SOL.lambdas],'-o')
        xlabel('Ux')
        ylabel('\lambda')
        grid on
        subplot(1,2,2)
        plot(U_plot_y, [0; SOL.lambdas],'-o')
        xlabel('Uy')
        ylabel('\lambda')
        grid on
    otherwise
        error('Solution type is not valid. Choose between: \n- Linear\n- Nonlinear')
end

sgtitle('Load path plots') 