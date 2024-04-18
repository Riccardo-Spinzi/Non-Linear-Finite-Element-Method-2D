function POST = POST_update_deformed_shape (POST, MODEL, nincr)

% --- Material coordinates stay the same at every load increment
POST.STEP(nincr).X = MODEL.XY(:,1);
POST.STEP(nincr).Y = MODEL.XY(:,2);

% --- Update spatial coordinates 
POST.STEP(nincr).x = MODEL.xy(:,1);
POST.STEP(nincr).y = MODEL.xy(:,2);

% --- Divide displacement components
U_res = reshape(MODEL.U,2,length(MODEL.U)/2)';
POST.STEP(nincr).Ux = U_res(:,1);
POST.STEP(nincr).Uy = U_res(:,2);

