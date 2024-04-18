function [ Bc, Bg ] = get_B_matrix( KINEMATICS, nl_check, formul_check )

switch lower(nl_check)
    case 'linear'
        Bc = [];
        Bg = [];

        dN=KINEMATICS.dN_spatial;

        for i = 1 : length(dN)
            Bc_temp = [dN(i,1)   0;
                  0         dN(i,2);
                  dN(i,2)   dN(i,1)];
            Bc = [Bc Bc_temp];
        end

    case 'nonlinear'

        Bc = [];
        Bg = [];
        

        % --- Check for Updated or total Lagrangian formulation
        if strcmpi(formul_check,'ul')
            
            dN=KINEMATICS.dN_spatial;
            for i = 1 : length(dN)
                Bc_temp = [dN(i,1)   0;
                      0         dN(i,2);
                      dN(i,2)   dN(i,1)];
                Bc = [Bc Bc_temp];
            end

        elseif strcmpi(formul_check,'tl')

            F = KINEMATICS.Def_gradient;
            dN=KINEMATICS.dN_material;
            
            for i = 1 : length(dN)
                Bc_temp = [F(1,1)*dN(i,1)                   F(2,1)*dN(i,1);
                           F(1,2)*dN(i,2)                   F(2,2)*dN(i,2);
                           F(1,2)*dN(i,1)+ F(1,1)*dN(i,2)   F(2,1)*dN(i,2)+ F(2,2)*dN(i,1)];
                Bc = [Bc Bc_temp];
            end

        end
        
        % Bg does not change between UL and TL
        Bg_temp = dN;
        Bg = [Bg Bg_temp']; 

    otherwise

        error('Unknown method')

end