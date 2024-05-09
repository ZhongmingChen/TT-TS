function [G, err, rel, T] = tt_als(Y, R, sigma, varargin) 

%% Include relevant files

addpath(genpath('help_functions'));

%% Handle optional inputs
params = inputParser;
addParameter(params, 'tol', 1e-12 , @isscalar);
addParameter(params, 'maxiters', 100, @(x) isscalar(x) & x > 0);
addParameter(params, 'verbose', false, @isscalar);
parse(params, varargin{:});
tol = params.Results.tol;
maxiters = params.Results.maxiters;

%% Initialize core tensors
    sizeY = size(Y); N = length(sizeY);
    G = cell(N,1); %rng(100);
    for n = 1:N
    %     G{n} = rand(R(n), R(n+1)*sizeY(n));    % R=[1 r r ... r 1];  
          G{n} = zeros(R(n) * R(n+1), sizeY(n));
    %     G{n} = eye(R(n), R(n+1)*sizeY(n));
    end

%% Main loop: Iterate until convergence, for a maximum of maxiters iterations
    err = []; 
    vecY = reshape(double(Y),[],1);
    rel = []; T = [];
    YsT = cell(N, 1); 
    Gs_hat = cell(N+1,1);  Gs_hat{1} = 1;  Gs_hat{N+1} = 1; 
    ite =1;    rel_err = tol +1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for n = 1:N
%         temp = double(tenmat(Y,n)).';
%         I = sort((randperm(prod(sizeY)/sizeY(n),J1)));
%         YsT{n} = temp(I,:);
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (ite <= maxiters)  &&  (rel_err > tol)
  
% generate Yst
    for n = 1:N
        YsT{n} = double(tenmat(Y,n)).';
    end    
    
% right-orthogonalization and generate Gs_hat
    for ind=N:-1:2
        temp1 = reshape(G{ind}, R(ind), R(ind+1)*sizeY(ind));
        [Q, Rfac] = qr(temp1.',0);        
        G{ind-1} = permute(reshape(G{ind-1}, R(ind-1), R(ind), sizeY(ind-1)), [1 3 2]); 
        G{ind-1} = reshape(G{ind-1}, R(ind-1)*sizeY(ind-1), R(ind))*Rfac.';
        R(ind) = size(Q, 2);    
        G{ind-1} = permute(reshape(G{ind-1}, R(ind-1), sizeY(ind-1), R(ind)), [1 3 2]);
        G{ind-1} = reshape(G{ind-1}, R(ind-1)*R(ind), sizeY(ind-1));
        G{ind} = Q.';    G{ind} = permute(reshape(G{ind}, R(ind),R(ind+1), sizeY(ind)),[1,3,2]);
        temp2 = reshape(G{ind},R(ind)*sizeY(ind),R(ind+1));
        G{ind} = reshape(permute(G{ind},[1,3,2]),[],sizeY(ind));  % Prepare for running tt_contraction(G)
        Gs_hat{ind} = temp2 * Gs_hat{ind+1};
        Gs_hat{ind} = reshape(Gs_hat{ind},R(ind),[]);
    end
    

    
% inner update
    err_core = 0;
    tic;
    for ind = 1:N
        H = kron(Gs_hat{ind+1}.',Gs_hat{ind});  B = YsT{ind};
        if sigma == 0
%            y = pinv((H'*H + sigma*eye(R(ind)*R(ind+1)))) * (H'*B);
%            y = (H'*H + sigma*eye(R(ind)*R(ind+1))) \ (H'*B);
             y = H' * B;
        else
            y = (H'*B+sigma*G{ind})/(1+sigma);
%           y = ((H'*H + sigma*eye(R(ind)*R(ind+1)))) \ (H'*B + sigma*G{ind});
        end
%         y = (SH'*SH + sigma*eye(R(ind)*R(ind+1))) \ (SH'*B);
%         y = (H'*H) \ (H'*B);
%         y = H \ B;
%         myfun = @(x) kron(eye(sizeY(ind)), (H'*H)) * x;    vecHtB = H'*B;    y = reshape(pcg(myfun, vecHtB(:)), size(G{ind}));

        err_core = max(err_core, norm(y - G{ind},'fro')/norm(y,'fro')); 

% Update Gs_hat to prepare the next iteration
        if ind < N
            yy = permute(reshape(y, R(ind), R(ind+1), sizeY(ind)), [1 3 2]);
            yy = reshape(yy, R(ind)*sizeY(ind), R(ind+1));
            [Q, Rfac] = qr(yy, 0);     
            G{ind+1} = Rfac*reshape(G{ind+1}, R(ind+1), R(ind+2)*sizeY(ind+1));
            R(ind+1) = size(Q, 2);    G{ind+1} = reshape(G{ind+1}, R(ind+1)*R(ind+2), sizeY(ind+1));
            G{ind} = reshape(Q, R(ind), sizeY(ind), R(ind+1));
            temp = reshape(G{ind},R(ind),[]);
            G{ind} = reshape(permute(G{ind},[1,3,2]),[],sizeY(ind)); %% Prepare for running tt_contraction(G)
            Gs_hat{ind+1} = Gs_hat{ind}*temp;
            Gs_hat{ind+1} = reshape(Gs_hat{ind+1},[],R(ind+1));
        else
            G{N} = y;
            t = toc; T = [T t];
            vecG = tt_contraction(G);
            rel_err = norm(vecG-vecY)/norm(vecY);
        end
    end
    
    fprintf(' Ite %2d: err_core = %7.1e , rel_err = %7.1e \n', ite, err_core, rel_err);
    err = [err err_core];
    rel = [rel rel_err];
%   log_rel = log10(rel);
    ite = ite + 1;
    
end

end
