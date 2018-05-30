using IterativeSolvers
include("usymqr.jl")
using Base.Test
using LinearMaps

srand(123)
n = 15

function unsymmetric_problem(T, n)
    A = rand(T, n, n)
    x = ones(T, n)
    b = A * x
    A = A + n * I
    A, x, b
end

srand(123)
n = 15


@testset "Unsymmetric Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A, x, b = unsymmetric_problem(T, n)
    tol = sqrt(eps(real(T)))
    x0 = rand(T, n)

    x1, hist1 = usymqr(A, b, maxiter = 10n, tol = tol, log = true)
    x2, hist2 = usymqr!(x0, A, b, maxiter = 10n, tol = tol, log = true)

    @test isa(hist1, ConvergenceHistory)
    @test norm(b - A * x1) / norm(b) ≤ tol
    @test hist1.isconverged
    @test norm(b - A * x2) / norm(b) ≤ tol
    @test x2 == x0
end

# @testset "SparseMatrixCSC{$T}" for T in (Float32, Float64, Complex64, Complex128)
#     A = let
#         B = sprand(n, n, 2 / n)
#     end
#
#     x = ones(T, n)
#     b = A * x
#     tol = sqrt(eps(real(T)))
#
#     x_approx, hist = minres(A, b, maxiter = 10n, tol = tol, log = true)
#
#     @test norm(b - A * x_approx) / norm(b) ≤ tol
#     @test hist.isconverged
# end
#
# end
