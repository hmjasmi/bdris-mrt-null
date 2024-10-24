% This function is used to implement WF at BS
% Inputs: - W: BS-BD-RIS Channel
%         - H: BD-RIS-Uk Channel
%         - Theta: BD-RIS Scattering matrix
%         - P_max: Max power budget at BS
%         - Approach: "1" for built-in function
%                     "2" for Bazzi method
%                     "3" for Bruno method
%                     "4" for G. Levin G. Levin (2024). Water Filling Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/3592-water-filling-algorithm), MATLAB Central File Exchange. Retrieved October 17, 2024.
% Output: Precoder matrix for WF
function P = func_Prec_WF(W,H,Theta,P_max,N0,Approach)
E = H*Theta*W;
K = size(E,1);
g = abs(diag(E).').^2;
if(Approach==1) % Built-in function
    P = waterfill(P_max,N0./g).';
    P = diag(sqrt(P));

elseif(Approach==2) % Bazzi
    num_points = 1e5;
    tol = min(1e-5,N0*1e-2);
    mu_axis = linspace(tol,5,num_points);
    Pn = max(1./mu_axis - N0./g.',0);
    dualFunc = sum(log(1 + 1./N0*Pn.*repmat(g.',[1,num_points]))) - mu_axis.*(sum(Pn) - P_max);
    [ming,idx] = find(dualFunc==min(dualFunc));

    mu = mu_axis(idx);

    P = max(1./mu - N0./g.',0);
    if (sum(P) ~= P_max)
        P = P * (P_max / sum(P));
    end
    P = diag(sqrt(P));

elseif(Approach==3) % Bruno
    rho = 1/N0;

    % Initialize power allocation vector
    P = zeros(K, 1);

    % Sort channels in decreasing order
    [g_sorted, idx] = sort(g.', 'descend');

    % Iterative power allocation
    for i = 1:K
        % Remaining channels
        remaining_channels = K - i + 1;
        % mu is the water-level
        % Compute the water-filling level mu for the remaining channels
        sum_lambda_inv = sum(1 ./ (rho * g_sorted(1:remaining_channels)));
        mu = (1 / remaining_channels) * (1 + sum_lambda_inv);

        % Calculate power for each channel and apply water-filling formula
        for k = 1:remaining_channels
            p_k = max(mu - 1 / (rho * g_sorted(k)), 0);
            P(idx(k)) = p_k;  % Store power in the original order
        end

        % Check if the last channel has non-negative power
        if P(idx(remaining_channels)) <= 0
            P(idx(remaining_channels)) = 0;  % Set negative power to zero
            break;  % Exit the loop if power constraint is fully utilized
        end
    end

    % Normalize power to satisfy the total power constraint
    if (sum(P) ~= P_max)
        P = P * (P_max / sum(P));
    end

    P = diag(sqrt(P));

elseif(Approach==4) % G. Levin
    % Define a tolerane for the water-level loop
    tol = min(1e-5,N0*1e-2);
    % lambda = 1/mu which is the water-level
    % Initialize lambda based on the best channel
    lambda = min(N0./g) + P_max/K;
    % Compute the used power
    ptot = sum(max(lambda - N0./g,0));
    % Modify the water-level to use the whole power budget
    % via AO
    while abs(P_max - ptot)>tol
        lambda = lambda + (P_max - ptot)/K;
        ptot = sum(max(lambda - N0./g,0));
    end
    P = max(lambda - N0./g,0);
    if (sum(P) ~= P_max)
        P = P * (P_max / sum(P));
    end
    P = diag(sqrt(P));
end

