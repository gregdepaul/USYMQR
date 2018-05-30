export usymqr_iterable, usymqr, usymqr!

import Base.LinAlg: BLAS.axpy!, givensAlgorithm
import Base: start, next, done

import IterativeSolvers.zerox
import IterativeSolvers.reserve!
import IterativeSolvers.setconv
import IterativeSolvers.shrink!
import IterativeSolvers.nextiter!


mutable struct USYMQRIterable{matT, solT, vecT <: DenseVector, smallVecT <: DenseVector, rotT <: Number, realT <: Real}
    A::matT
    x::solT

    # Diagonalization Constants
    beta::realT
    gamma::realT
    gamma_prev::realT
    proj::realT

    # Krylov basis vectors
    v_prev::vecT
    v::vecT
    u::vecT

    # Plane Rotation Update Vectors
    p::vecT
    q::vecT

    # W = R * inv(V) is computed using 3-term recurrence
    w1::vecT
    w2::vecT

    # rhs is just two active values of the right-hand side.
    rhs1::smallVecT
    rhs2::smallVecT

    # Some Givens rotations
    c_prev::rotT
    c::rotT
    s::rotT

    # QR Constants
    tau::realT
    sbar::realT
    tau_prev::realT

    # Convergence Criterion
    q1::realT
    q2::realT
    AAnorm::realT

    # Bookkeeping
    mv_products::Int
    maxiter::Int
    tolerance::realT
    resnorm::realT
end

function usymqr_iterable!(x, A, b;
    initially_zero::Bool = false,
    tol = sqrt(eps(real(eltype(b)))),
    maxiter = size(A, 2)
)

    (m,n) = size(A);
    T = eltype(x)
    HessenbergT = real(T)

    c = rand(T, n, 1);

    beta1 = norm(b);
    gamma = norm(c);

    v_prev = similar(b)
    v = similar(b)
    copy!(v, b / beta1)

    u = similar(c)
    copy!(u, c / gamma)

    p = similar(v);
    q = similar(u);
    beta = 0;
    gamma = 0;

    sbar = gamma;
    rhs1 = beta1;
    AAnorm = 0;
    tau = 0;
    gamma_prev = 0;
    c_prev = 1;
    c = 1;
    s = 0;
    q1 = 0;
    q2 = 0;
    proj = 0;

    w1 = similar(x)
    w2 = similar(x)

    mv_products = 0

    # For nonzero x's, we must do an MV for the initial residual vec
    if !initially_zero
        # Use v_next to store Ax; v_next will soon be overwritten.
        A_mul_B!(v_next, A, x)
        axpy!(-one(T), v_next, v_curr)
        mv_products = 1
    end

    resnorm = beta1;
    reltol = resnorm * tol;

    USYMQRIterable(
        A, x,
        beta, gamma, gamma_prev, proj,
        v_prev, v, u,
        p, q,
        w1, w2,
        rhs1, 0.0,
        c_prev, c, s,
        tau, sbar, 0.0,
        q1, q2, AAnorm,
        mv_products, maxiter, reltol, resnorm
    )
end

Arnorm(m::USYMQRIterable) = m.resnorm*norm([m.gamma_prev*m.q1 + m.proj*m.q2; m.gamma*m.q2]);
Anorm(m::USYMQRIterable) = sqrt(m.AAnorm);
AAnorm(m::USYMQRIterable) = AAnorm(m) + m.proj^2 + m.beta^2 + m.gamma^2;


condition1(m::USYMQRIterable) = (Arnorm(m)/(Anorm(m)*m.resnorm) < m.tolerance);

condition2(m::USYMQRIterable) = (m.resnorm < m.tolerance*Anorm(m) + m.tolerance);

converged(m::USYMQRIterable) = m.resnorm ≤ m.tolerance

start(::USYMQRIterable) = 1

done(m::USYMQRIterable, iteration::Int) = iteration > m.maxiter || condition1(m) || condition2(m);

function next(m::USYMQRIterable, iteration::Int)

    m.v_prev = m.v;
    m.tau_prev = m.tau;

    # p = A*v - gama*p;
    p = A_mul_B(p, m.A, m.v);
    axpy!(-m.gamma, m.p, p);
    m.p = p;

    #q = A'*u - beta*q;
    q = Ac_mul_B(q, m.A, m.u);
    axpy!(-m.beta, m.q, q);
    m.q = q;

    # Orthogonalize w.r.t. m.u
    m.proj = dot(m.p, m.u);
    temp = m.u;
    m.u = m.p;
    axpy!(-m.proj, temp, m.u)
    m.p = temp;

    # Orthogonalize w.r.t. m.v
    temp = m.v;
    m.v = m.q;
    axpy!(-m.proj, temp, m.v)
    m.q = temp;

    # Normalize u, v with beta and gamma
    m.beta = norm(m.u);
    m.gamma = norm(m.v);
    m.u /= m.beta;
    m.v /= m.gamma;

    # Consider using c, s, m.H[3] = givensAlgorithm(m.H[3], m.H[4])
    # Form QR factorization
    sigma = m.c * m.sbar + m.s * proj;
    rbar = -m.s * m.sbar + m.c * proj;
    m.tau = m.s * m.gamma;
    m.sbar = m.c * m.gama;

    if iteration == 1
        rbar = temp;
        m.sbar = m.gamma;
    end

    rho = sqrt(m.rbar^2 + m.beta^2);
    m.c = m.rbar/rho;
    m.s = m.beta/rho;

    # Next rotation
    temp = m.rhs1;
    m.rhs1 =  m.c * temp;
    m.rhs2 = -m.s * temp;

    # Update solution x
    updateSoln!(m.x, m.rhs1, m.v_prev, m.tau_prev, sigma, rho, m.w1, m.w2)

    # The approximate residual is cheaply available
    m.resnorm = abs(m.rhs2);
    m.q1 = -m.c_prev * m.s;
    m.q2 = m.c;

    m.rhs1 = m.rhs2;
    m.gamma_prev = m.gamma;
    m.c_prev = m.c;

    m.resnorm, iteration + 1
end


function updateSoln!(x, rh_s1, v_prev, tau_prev, sigma, rho, w1, w2)
    # Update solution
    w3 = ( v_prev - sigma*w2 - tau_prev*w1 )/rho;
    x = x + rh_s1*w3;
    w1 = w2;
    w2 = w3;
end


"""
    usymqr!(x, A, b; kwargs...) -> x, [history]

Solve Ax = b for Unsymmetric matrices A using usymqr.

# Arguments

- `x`: initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `initially_zero::Bool = false`: if `true` assumes that `iszero(x)` so that one
  matrix-vector product can be saved when computing the initial
  residual vector;
- `tol`: tolerance for stopping condition `|r_k| / |r_0| ≤ tol`. Note that the residual is computed only approximately;
- `maxiter::Int = size(A, 2)`: maximum number of iterations;
- `log::Bool = false`: keep track of the residual norm in each iteration;
- `verbose::Bool = false`: print convergence information during the iterations.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;
- `history`: convergence history.
"""
function usymqr!(x, A, b;
    verbose::Bool = false,
    log::Bool = false,
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Int = size(A, 2),
    initially_zero::Bool = false
)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && reserve!(history, :resnorm, maxiter)

    iterable = usymqr_iterable!(x, A, b;
        tol = tol,
        maxiter = maxiter,
        initially_zero = initially_zero
    )

    if log
        history.mvps = iterable.mv_products
    end

    # for (iteration, resnorm) = enumerate(iterable)
    #     if log
    #         nextiter!(history, mvps = 1)
    #         push!(history, :resnorm, resnorm)
    #     end
    #     verbose && @printf("%3d\t%1.2e\n", iteration, resnorm)
    # end

    verbose && println()
    log && setconv(history, converged(iterable))
    log && shrink!(history)

    log ? (iterable.x, history) : iterable.x
end

"""
    usymqr(A, b; kwargs...) -> x, [history]

Same as [`usymqr!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
usymqr(A, b; kwargs...) = usymqr!(zerox(A, b), A, b; initially_zero = true, kwargs...)
