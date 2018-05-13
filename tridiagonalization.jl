function tridiagonalization(A, epsilon, max_iterations)

    (m,m) = size(A);
    b = randn(1 ,m);
    c = randn(1, m);

    beta = [];
    gamma = [];
    alpha = [];

    push!(beta, norm(b));
    push!(gamma, norm(c));

    P = zeros(m,max_iterations);
    Q = zeros(m,max_iterations);

    P[:, 1] = b / beta[1]
    Q[:, 1] = c / gamma[1]

    i = 1;
    while i < max_iterations

        u = A*Q[:, i];
        v = A'*P[:, i];

        if i > 1
            u = u - gamma[i]*P[:, i - 1];
            v = v - beta[i]*Q[:, i - 1];
        end

        push!(alpha, P[:, i]'*u);

        u = u - alpha[i]*P[: ,i];
        v = v - alpha[i]*Q[:, i];

        if (norm(u) < epsilon || norm(v) < epsilon)
            break;
        else
            push!(beta, norm(u));
            push!(gamma, norm(v));

            P[:, i+1] = u/beta[i+1];
            Q[:, i+1] = v/gamma[i+1];
        end

        i = i + 1;
    end

    betaFT = convert(Array{Float64,1}, beta);
    alphaFT = convert(Array{Float64,1}, alpha);
    gammaFT = convert(Array{Float64,1}, gamma);

    T = Tridiagonal(betaFT[2:size(alphaFT,1)], alphaFT, gammaFT[2:size(alphaFT,1)]);

    P = P[:, 1:length(alphaFT)];
    Q = Q[:, 1:length(alphaFT)];

    return P, T, Q;

    end
