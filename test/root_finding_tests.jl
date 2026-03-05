using Test
using TurboMachineModel
using Random

@testset "Root Finding" begin
    U = TurboMachineModel.Utility

    @testset "Polynomial Roots" begin
        # (x-1)*(x-2)*(x-3) = x^3 - 6x^2 + 11x - 6
        f_poly(x) = x^3 - 6x^2 + 11x - 6
        roots = U.bracket_bisect_roots(f_poly, (0.5, 3.5); n_scan=301, root_tol=1e-10)
        @test length(roots) == 3
        @test isapprox(roots[1], 1.0; atol=1e-8)
        @test isapprox(roots[2], 2.0; atol=1e-8)
        @test isapprox(roots[3], 3.0; atol=1e-8)
    end

    @testset "Exponential Root" begin
        f_exp(x) = exp(x) - 2.0
        roots = U.bracket_bisect_roots(f_exp, (-2.0, 2.0); n_scan=201, root_tol=1e-10)
        @test length(roots) == 1
        @test isapprox(roots[1], log(2.0); atol=1e-8)
    end

    @testset "Rational With Non-Finite Region" begin
        # Asymptote at x = 1.0; finite roots at x = 0 and x = 2
        f_rat(x) = (x * (x - 2.0)) / (x - 1.0)
        roots = U.bracket_bisect_roots(f_rat, (-0.5, 2.5); n_scan=401, root_tol=1e-10)
        @test length(roots) == 2
        @test isapprox(roots[1], 0.0; atol=1e-8)
        @test isapprox(roots[2], 2.0; atol=1e-8)
    end

    @testset "No Root / Invalid Domain Handling" begin
        f_no_root(x) = x^2 + 1.0
        roots = U.bracket_bisect_roots(f_no_root, (-2.0, 2.0))
        @test isempty(roots)

        f_nan(x) = x > 0 ? NaN : x + 1.0
        roots_nan = U.bracket_bisect_roots(f_nan, (-2.0, 2.0))
        @test length(roots_nan) == 1
        @test isapprox(roots_nan[1], -1.0; atol=1e-8)
    end

    @testset "Continuation Hints and Dedupe" begin
        f_poly(x) = x^3 - 6x^2 + 11x - 6
        roots = U.bracket_bisect_roots(
            f_poly,
            (0.5, 3.5);
            n_scan=31, # intentionally coarse
            prior_roots=[1.0, 2.0, 3.0],
            root_tol=1e-10,
        )
        @test length(roots) == 3
        @test isapprox(roots[1], 1.0; atol=1e-8)
        @test isapprox(roots[2], 2.0; atol=1e-8)
        @test isapprox(roots[3], 3.0; atol=1e-8)

        # Root very close to a probe point can produce near-duplicates.
        f_near(x) = (x - 0.5) * (x - 0.50000001)
        roots_default = U.bracket_bisect_roots(f_near, (0.0, 1.0); n_scan=401, root_tol=1e-12)
        roots_merged = U.bracket_bisect_roots(
            f_near,
            (0.0, 1.0);
            n_scan=401,
            root_tol=1e-12,
            dedupe_atol=1e-5,
        )
        @test length(roots_default) >= 1
        @test length(roots_merged) == 1
        @test isapprox(roots_merged[1], 0.5; atol=1e-4)
    end

    @testset "Input Validation" begin
        @test_throws ErrorException U.bracket_bisect_roots(x -> x, (1.0, 1.0))
        @test_throws ErrorException U.bracket_bisect_roots(x -> x, (0.0, 1.0); n_scan=3)
    end

    @testset "Feasibility Backoff" begin
        eval_fn(v) = (value=v, feasible=(v <= 4.0))
        is_feasible(r) = r.feasible

        hit = U.feasibility_backoff(
            eval_fn,
            3.5;
            enabled=true,
            min_value=0.0,
            max_value=5.0,
            is_feasible=is_feasible,
        )
        @test hit.converged
        @test !hit.used_backoff
        @test isapprox(hit.value, 3.5; atol=1e-12)

        backed = U.feasibility_backoff(
            eval_fn,
            5.0;
            enabled=true,
            min_value=0.0,
            max_value=5.0,
            value_tol=1e-6,
            max_iters=40,
            n_probe=41,
            is_feasible=is_feasible,
        )
        @test backed.converged
        @test backed.used_backoff
        @test isapprox(backed.value, 4.0; atol=2e-2)
        @test backed.result.feasible

        disabled = U.feasibility_backoff(
            eval_fn,
            5.0;
            enabled=false,
            min_value=0.0,
            max_value=5.0,
            is_feasible=is_feasible,
        )
        @test !disabled.converged
        @test isnan(disabled.value)
        @test disabled.result === nothing

        none = U.feasibility_backoff(
            v -> (value=v, feasible=false),
            5.0;
            enabled=true,
            min_value=0.0,
            max_value=5.0,
            is_feasible=is_feasible,
        )
        @test !none.converged
        @test isnan(none.value)
        @test none.result === nothing
    end

    @testset "Randomized Polynomial Roots" begin
        rng = MersenneTwister(20260305)

        # Property-style test: for randomly generated cubic polynomials with
        # known distinct real roots in the scan interval, the root finder should
        # recover all roots within tolerance.
        for _ in 1:50
            roots_true = sort(rand(rng, 3) .* 5 .- 2.5) # roots in [-2.5, 2.5]
            # enforce minimum spacing to avoid near-coincident roots collapsing
            if minimum(diff(roots_true)) < 0.12
                continue
            end

            scale = 10.0^(rand(rng) * 4 - 2) # scale in [1e-2, 1e2]
            f(x) = scale * (x - roots_true[1]) * (x - roots_true[2]) * (x - roots_true[3])

            roots_found = U.bracket_bisect_roots(
                f,
                (-3.0, 3.0);
                n_scan=801,
                root_tol=1e-10,
                max_bisect_iters=120,
                dedupe_atol=1e-7,
            )

            @test length(roots_found) == 3
            @test isapprox(roots_found[1], roots_true[1]; atol=5e-6)
            @test isapprox(roots_found[2], roots_true[2]; atol=5e-6)
            @test isapprox(roots_found[3], roots_true[3]; atol=5e-6)
        end
    end
end
