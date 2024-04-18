function POST = POST_strain_recovery( MODEL, POST, KINEMATICS, nincr )

% --- Save deformation gradients
POST.STEP(nincr).F= KINEMATICS.F;

[~,l] = size(MODEL.nodes);
nint_p = length(MODEL.int_rule.x);

for i = 1 : MODEL.nels

    for j = 1 : nint_p^l

        F = POST.STEP(nincr).F(i,j).F_mat;
        C = F' * F;
        I = eye( l, l );

        % --- Principal stretches and RU/VR decomposition
        [ N, Lambda2 ] = eig( C );
        Lambda = sqrt( Lambda2 );
        U = N * Lambda * N';
        R = U \ F;
        V = F * R';
        n = R * N;

        % --- Green Lagrange & Almansi strains
        POST.STEP(nincr).E2(i,j).E = 0.5 * N * ( Lambda2 - I ) * N';
        Lambdam2 = diag( 1 ./ diag( Lambda2 ) );
        POST.STEP(nincr).e2(i,j).e = 0.5 * n * ( I - Lambdam2 ) * n' ;

        % --- Logarithmic strains
        logLambda = diag( log( diag( Lambda ) ) );
        POST.STEP(nincr).E0(i,j).E = N * logLambda * N';
        POST.STEP(nincr).e0(i,j).e = n * logLambda * n'; 

    end
end

