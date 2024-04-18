function MODEL = stress_plot( MODEL, POST, SOL )

Barlow = MODEL.Barlow_points;
nl_check = SOL.type;

if strcmpi(nl_check,'linear') == 1 
    [m,n] = size(MODEL.elements);
    X = zeros(m,n);
    Y  = zeros(m,n);
    
    for i = 1 : MODEL.nels
        X(i,:) = MODEL.nodes(MODEL.elements(i,:),1);
        Y(i,:) = MODEL.nodes(MODEL.elements(i,:),2);
    end

    if (n == 8 && Barlow == 0) || n == 8
        
        X(:,[5:end]) = [];
        Y(:,[5:end]) = [];
    
    end

    figure
    colormap("turbo")
    subplot(1,3,1)
    patch(X',Y',MODEL.sigma_xx')
    xlabel('x')
    ylabel('y')
    title('\sigma_{xx}')
    colorbar
    subplot(1,3,2)
    patch(X',Y',MODEL.sigma_yy')
    xlabel('x')
    ylabel('y')
    title('\sigma_{yy}')
    colorbar
    subplot(1,3,3)
    patch(X',Y',MODEL.sigma_xy')
    xlabel('x')
    ylabel('y')
    title('\sigma_{xy}')
    colorbar

else
    % --- NONLINEAR STRESS PLOT
    [m,n] = size(MODEL.elements);
    sigma_xx = zeros(MODEL.nels,length(MODEL.int_rule.x)^2);
    sigma_yy = zeros(size(sigma_xx));
    sigma_xy = zeros(size(sigma_xx));

    if n > 4
        MODEL.elements(:,5:end) = [];
        n = 4;
    end
    
    x = zeros(m,n);
    y = zeros(m,n);

    for i = 1 : MODEL.nels

        cont = 1;

        for j = 1 : length(MODEL.int_rule.x)^2

            if (j == 1 || j == 3 || j == 7 || j==9) && n > 4

                sigma_xx(i,cont) = POST.STEP(end).stress(i,j).sigma(1,1);
                sigma_yy(i,cont) = POST.STEP(end).stress(i,j).sigma(2,2);
                sigma_xy(i,cont) = POST.STEP(end).stress(i,j).sigma(1,2);

                cont = cont + 1;

            else

                sigma_xx(i,j) = POST.STEP(end).stress(i,j).sigma(1,1);
                sigma_yy(i,j) = POST.STEP(end).stress(i,j).sigma(2,2);
                sigma_xy(i,j) = POST.STEP(end).stress(i,j).sigma(1,2);

            end

        end

        x(i,:) = POST.STEP(end).x(MODEL.elements(i,:));
        y(i,:) = POST.STEP(end).y(MODEL.elements(i,:));

    end
    
    sigma_xx(:,5:end) = [];
    sigma_yy(:,5:end) = [];
    sigma_xy(:,5:end) = [];

     % MANIPULATE STRESSES FOR PATCH command
%     sigma_xx(:,[2,3]) = sigma_xx(:,[3,2]);
    sigma_xx(:,[3,4]) = sigma_xx(:,[4,3]);
%     sigma_yy(:,[2,3]) = sigma_yy(:,[3,2]);
    sigma_yy(:,[3,4]) = sigma_yy(:,[4,3]);
%     sigma_xy(:,[2,3]) = sigma_xy(:,[3,2]);
    sigma_xy(:,[3,4]) = sigma_xy(:,[4,3]);

    MODEL.sigma_xx = sigma_xx;
    MODEL.sigma_yy = sigma_yy;
    MODEL.sigma_xy = sigma_xy;

    figure(2)
    colormap("turbo")
    subplot(1,3,1)
    patch(x',y',sigma_xx')
    xlabel('x')
    ylabel('y')
    title('\sigma_{xx}')
    colorbar
    subplot(1,3,2)
    patch(x',y',sigma_yy')
    xlabel('x')
    ylabel('y')
    title('\sigma_{yy}')  
    colorbar
    subplot(1,3,3)
    patch(x',y',sigma_xy')
    xlabel('x')
    ylabel('y')
    title('\sigma_{xy}')
    colorbar
end
