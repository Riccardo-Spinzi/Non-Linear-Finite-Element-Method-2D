function POST = POST_update_deformed_shape (POST, MODEL, nincr)

         % % --------------- FUNCTION INFO ---------------- % %

% POST_update_deformed_shape updates the shape of the structure at the end
% of every load step. 
%
%         POST = POST_update_deformed_shape (POST, MODEL, nincr)
%
% -------------------------------------------------------------------------
% Input arguments: 
% POST                [struct]      POST structure                  [multi]
% MODEL               [struct]      MODEL structure                 [multi]
% nincr               [1x1 double]  load step number                [-]            
%
% -------------------------------------------------------------------------
% Output arguments:
% POST                [struct]      struct containing parameters     
%                                   for the POST of the structure   [multi]
%
% -------------------------------------------------------------------------

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

