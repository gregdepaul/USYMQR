export usymqr_iterable, usymqr, usymqr!

import Base.LinAlg: BLAS.axpy!, givensAlgorithm
import Base: start, next, done

mutable struct USYMQRIterable{matT, solT, vecT <: DenseVector, smallVecT <: DenseVector, rotT <: Number, realT <: Real}
    A::matT
    x::solT

    beta1::realT
    gama1::realT
    beta::realT
    gama::realT

    AAnorm::realT
    sbar::realT
    rhs1::realT
    tau::realT
    gamold::realT
    cold::realT

    # Krylov basis vectors
    u::vecT
    v::vecT
    p::vecT
    q::vecT

    # W = R * inv(V) is computed using 3-term recurrence
    w1::vecT
    w2::vecT
    q1::vecT
    q2::vecT

    # Vector of size 4, holding the active column of the Hessenberg matrix
    # rhs is just two active values of the right-hand side.
    H::smallVecT
    rhs::smallVecT

    # Some Givens rotations
    c::rotT
    s::rotT

    # Bookkeeping
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
    c = randn(1, m);

    beta1  = norm(b);
    u     = b/beta1;
    gama1 = norm(c);
    v     = c/gama1;

    p = zeros(m,1);
    q = zeros(n,1);
    beta = 0;
    gama = 0;

    AAnorm = 0;
    sbar   = gama;
    rhs1   = beta1;
    tau    = 0;
    gamold = 0;
    cold = 1;
    c = 1; s = 0;
    resnorm = beta1;
    q1 = 0; q2 = 0;

    w1 = zeros(n,1);
    w2 = zeros(n,1);

    reltol = resnorm * tol;

    USYMQRIterable(
        A, x,
        beta1, u, gama1, v,
        p, q, beta, gama,
        AAnorm, sbar, rhs1,tau, gamold, cold,
        c, s,
        q1, q2,
        w1, w2,
        maxiter, reltol, resnorm
    )
end

converged(m::USYMQRIterable) = m.resnorm ≤ m.tolerance

start(::USYMQRIterable) = 1

done(m::USYMQRIterable, iteration::Int) = iteration > m.maxiter || converged(m)

function next(m::USYMQRIterable, iteration::Int)
    # v_next = A * v_curr - H[2] * v_prev
    A_mul_B!(m.p, m.A, m.v_curr)

    iteration > 1 && axpy!(-m.H[2], m.v_prev, m.v_next)

    # Orthogonalize w.r.t. v_curr
    proj = dot(m.v_curr, m.v_next)
    m.H[3] = real(proj)
    axpy!(-proj, m.v_curr, m.v_next)

    # Normalize
    m.H[4] = norm(m.v_next)
    scale!(m.v_next, inv(m.H[4]))

    # Rotation on H[1] and H[2]. Note that H[1] = 0 initially
    if iteration > 2
        m.H[1] = m.s_prev * m.H[2]
        m.H[2] = m.c_prev * m.H[2]
    end

    # Rotation on H[2] and H[3]
    if iteration > 1
        tmp = -conj(m.s_curr) * m.H[2] + m.c_curr * m.H[3]
        m.H[2] = m.c_curr * m.H[2] + m.s_curr * m.H[3]
        m.H[3] = tmp
    end

    # Next rotation
    c, s, m.H[3] = givensAlgorithm(m.H[3], m.H[4])

    # Apply as well to the right-hand side
    m.rhs[2] = -conj(s) * m.rhs[1]
    m.rhs[1] = c * m.rhs[1]

    # Update W = V * inv(R). Two axpy's can maybe be one MV.
    copy!(m.w_next, m.v_curr)
    iteration > 1 && axpy!(-m.H[2], m.w_curr, m.w_next)
    iteration > 2 && axpy!(-m.H[1], m.w_prev, m.w_next)
    scale!(m.w_next, inv(m.H[3]))

    # Update solution x
    axpy!(m.rhs[1], m.w_next, m.x)

    # Move on: next -> curr, curr -> prev
    m.v_prev, m.v_curr, m.v_next = m.v_curr, m.v_next, m.v_prev
    m.w_prev, m.w_curr, m.w_next = m.w_curr, m.w_next, m.w_prev
    m.c_prev, m.s_prev, m.c_curr, m.s_curr = m.c_curr, m.s_curr, c, s
    m.rhs[1] = m.rhs[2]

    # Due to symmetry of the tri-diagonal matrix
    m.H[2] = m.H[4]

    # The approximate residual is cheaply available
    m.resnorm = abs(m.rhs[2])

    m.resnorm, iteration + 1
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

    for (iteration, resnorm) = enumerate(iterable)
        if log
            nextiter!(history, mvps = 1)
            push!(history, :resnorm, resnorm)
        end
        verbose && @printf("%3d\t%1.2e\n", iteration, resnorm)
    end

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
