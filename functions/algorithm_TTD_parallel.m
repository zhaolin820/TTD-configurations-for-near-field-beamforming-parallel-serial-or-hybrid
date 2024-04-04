function [R_convergence] = algorithm_TTD_parallel(para, h, P, user_r, user_theta)
%Penalty-based algorithm for parallel configurations
%  [R_convergence] = algorithm_TTD_parallel(para, h, P, user_r, user_theta)
%Inputs:
%   para: structure of the initial parameters
%   h: channel for all users
%   P: initial beamforming vectors
%   user_r: distance for all users (for initialization)
%   user_theta: angle for all users (for initialization)
%Outputs:
%   R_convergence: achievable rates at each iteration
%Date: 04/04/2024


%% initialization
t = zeros(para.N_T, para.N_RF);
A_PS = zeros(para.N, para.N_RF);
for n = 1:para.N_RF
    [a_n, t_n] = initialize_analog(para, user_theta(n), user_r(n));
    t(:,n) = t_n;
    A_PS(:,n) = a_n;
end

A = zeros(para.N, para.N_RF, para.M);
D = zeros(para.N_RF, para.K, para.M);
for m = 1:para.M
    A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m));
    D(:,:,m) = pinv(A(:,:,m))*P(:,:,m);
end

% R_convergence = zeros(1,100);
i = 1;
%% Optimization
rho = 1e4;
for outer_step = 1:100
    obj_pre = 0;
    for inner_step = 1:40
        [P] = update_auxiliary_var(para, h, P, A, D, rho);
        [A, t] = update_analog(para, P, D, t);
        [D] = update_digital(para, P, A, D);
        
        % objective value
        [R_sum_FD] = rate_fully_digital(para, P, h);
        penalty = 0;
        for m = 1:para.M
            Pm = P(:,:,m); Am = A(:,:,m); Dm = D(:,:,m);
            penalty = penalty + norm(Pm*pinv(Dm) - Am, 'fro')^2;
        end
        % check convergence of inner loops
        obj = R_sum_FD - 1/rho*penalty;
        if abs((obj-obj_pre)/obj) < 1e-3
            break;
        end
        obj_pre = obj;

        Pa = zeros(para.N, para.K, para.M);
        for m = 1:para.M
            Pa(:,:,m) = A(:,:,m)*D(:,:,m);
        end
        [R_sum] = rate_fully_digital(para, Pa, h);

        R_convergence(i) = R_sum/(para.M+para.Lcp); i = i+1;

        disp(['inner loop - ' num2str(inner_step) ', obj - ' num2str(obj) ', rate - ' num2str(R_sum_FD/(para.M+para.Lcp))...
            ', rate_HB - ' num2str(R_sum/(para.M+para.Lcp)) ', penalty - ' num2str(penalty)]);
    end
    [xi] = violation(para, P, A, D);    

    rho = 0.1*rho;
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% outer loop - ' num2str(outer_step) ', violation - ' num2str(xi)...
        ', rho - ' num2str(rho) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    if xi < 1e-5
        break;
    end
end

end

%% Update auxiliary variables
function [P] = update_auxiliary_var(para, h, P, A, D, rho)
    E = eye(para.K);
    for m = 1:para.M
        hm = h(:,:,m); Dm = D(:,:,m); Am = A(:,:,m); Pm = P(:,:,m);
        Phi = 0; Upsilon = 0;
        for k = 1:para.K
            hk = hm(:,k);
            pk = Pm(:,k); 
            I = norm(hk'*Pm)^2 + norm(Pm, 'fro')^2/para.Pt; 
            w_k = 1 + abs(hk'*pk)^2 / (I - abs(hk'*pk)^2);
            v_k = hk'*pk / I;
            
            Phi = Phi + w_k*abs(v_k)^2 * ( hk*hk' + eye(para.N)/para.Pt );
            Upsilon = Upsilon + w_k*conj(v_k)*E(:,k)*hk';
        end
        Upsilon = Upsilon + 1/rho * pinv(Dm)*Am';
        Psi = 1/rho * pinv(Dm)*pinv(Dm)';
        Pm = sylvester(Psi, Phi, Upsilon);
        P(:,:,m) = Pm';
    end
end

%% Update analog beamformer
function [A, t] = update_analog(para, P, D, t)

    t_search = 0:para.t_max/1e3:para.t_max;
    PD = zeros(para.N, para.N_RF, para.M);
    for m = 1:para.M
        PD(:,:,m) = P(:,:,m)*pinv(D(:,:,m));
    end

    % update PSs
    N_sub = para.N/para.N_T;
    A_PS = zeros(para.N, para.N_RF);
    for n = 1:para.N_RF
        for q = 1:para.N_T
            p_nq = PD((q-1)*N_sub+1:q*N_sub, n, :);  
            p_nq = squeeze(p_nq);
            t_nq = t(q,n);
            a_nq = sum(p_nq*diag(exp(1i*2*pi*para.fm_all*t_nq)), 2);
            A_PS((q-1)*N_sub+1:q*N_sub, n) = a_nq./abs(a_nq);
        end
    end

    % update TTDs
    t = zeros(para.N_T, para.N_RF);
    for n = 1:para.N_RF
        for q = 1:para.N_T
            p_nq = PD((q-1)*N_sub+1:q*N_sub, n, :);  
            p_nq = squeeze(p_nq);
            a_nq = A_PS((q-1)*N_sub+1:q*N_sub, n);
            psi_nq = p_nq'*a_nq;
            results = real(psi_nq.'*exp(-1i*2*pi*para.fm_all'*t_search));
            [~,I] = max(results);
            t(q,n) = t_search(I);
        end
    end

    % calculate overall analog beamformer
    A = zeros(para.N, para.N_RF, para.M);
    for m = 1:para.M
        A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m));
    end
end

%% Update digital beamformer
function [D] = update_digital(para, P, A, D)
    for m = 1:para.M
        D(:,:,m) = pinv(A(:,:,m))*P(:,:,m);
    end
end


%% calculate the overall analog beamformer
function [Am] = analog_bamformer(para, A, t, fm)
    e = ones(para.N/para.N_T,1);
    Tm = exp(-1i*2*pi*fm*t);
    Am = A .* kron(Tm, e);
end


%% constraint violation
function [xi] = violation(para, P, A, D)
    xi = zeros(para.M,1);
    for m = 1:para.M
        Pm = P(:,:,m); Am = A(:,:,m); Dm = D(:,:,m);
        xi(m) = norm(Pm - Am*Dm, 'inf');
    end
    xi = max(xi);
end


%% initialization analog beamformer
function [a, t] = initialize_analog(para, theta, r)
    c = 3e8;
    N_sub = para.N/para.N_T;
    t = zeros(para.N_T, 1);
    a = zeros(para.N, 1);
    for q = 1:para.N_T
        xi_q = (q-1-(para.N_T-1)/2)*N_sub;
        r_q = sqrt(r^2 + xi_q^2*para.d^2 - 2*r*xi_q*para.d*cos(theta));
        theta_q = acos( (r*cos(theta) - xi_q*para.d)/r_q );
        [beta, r_ave] = beta_func(r_q, theta_q, N_sub, para.d, para.fc);
        % time delay
        t(q) = - (r_q - r + r_ave)/c;
    
        % PS beamformer
        a((q-1)*N_sub+1 : q*N_sub) = conj(beta);
    end
    t = t - min(t);
    t(t>para.t_max) = para.t_max;
end
