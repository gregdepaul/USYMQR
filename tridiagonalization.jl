A = rand(4, 4)

(m,n) = size(A);
b = randn(1 ,n);
c = randn(1, m);

beta = [];
gamma = [];
alpha = [];

push!(beta, norm(b));
push!(gamma, norm(c));

P = zeros(m,n);
Q = zeros(m,n);

P[:, 1] = b / beta[1]
Q[:, 1] = c / gamma[1]

i = 1;
u = A*Q[:, i];
v = ctranspose(A)*P[:, i];

maxiters = 5;
while i < maxiters
    if i > 1
        u = u - gamma[i]*P[:, i - 1];
        v = v - beta[i]*Q[:, i - 1];
    end

    push!(alpha, ctranspose(P[:, i])*u);

    #print(alpha)
    u = u - alpha[i]*P[: ,i];
    v = v - ctranspose(alpha[i])*Q[:, i];

    if (norm(u) == 0 || norm(v) == 0)
        break;
    else
        push!(beta, norm(u));
        push!(gamma, norm(v));

        P[:, i] = u/beta[i+1];
        Q[:, i] = v/gamma[i+1];
    end

    i = i + 1;
end



betaFT = convert(Array{Float64,1}, beta);
alphaFT = convert(Array{Float64,1}, alpha);
gammaFT = convert(Array{Float64,1}, gamma);

T = Tridiagonal(betaFT[2:size(betaFT,1)-1], alphaFT, gammaFT[2:size(gammaFT,1)-1]);
