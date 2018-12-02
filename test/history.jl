using IterativeSolvers
using RecipesBase
using Test

@testset "ConvergenceHistory" begin

    KW = Dict{Symbol, Any}

    RecipesBase.is_key_supported(k::Symbol) = k == :sep ? false : true

    begin
        history = ConvergenceHistory(partial = false)
        history.iters = 3
        history.isconverged = true
        @test string(history) == "Converged after 3 iterations."
        history.isconverged = false
        @test string(history) == "Not converged after 3 iterations."
    end

    # No plottables
    begin
        history = ConvergenceHistory(partial = false)
        @test_throws ErrorException RecipesBase.apply_recipe(KW(), history)
    end

    # Without restart markers
    begin
        history = ConvergenceHistory(partial = false)
        history.iters = 3
        history.data[:resnorm] = [10.0, 3.0, 0.1]

        plots = (
            RecipesBase.apply_recipe(KW(), history),
            RecipesBase.apply_recipe(KW(), history, :resnorm)
        )

        for data in plots
            @test length(data) == 1
            @test data[1].plotattributes[:label] == "resnorm"
            @test data[1].plotattributes[:seriestype] == :line
        end
    end

    # With restart markers
    begin
        history = ConvergenceHistory(partial = false, restart = 2)
        history.iters = 3
        history.data[:resnorm] = [10.0, 3.0, 0.1]

        plots = (
            RecipesBase.apply_recipe(KW(), history),
            RecipesBase.apply_recipe(KW(), history, :resnorm)
        )

        for data in plots
            @test length(data) == 2
            @test data[2].plotattributes[:linecolor] == :white
        end
    end

    # Custom color
    begin
        history = ConvergenceHistory(partial = false, restart = 2)
        history.iters = 3
        history.data[:resnorm] = [10.0, 3.0, 0.1]

        plots = (
            RecipesBase.apply_recipe(KW(:sep => :red), history),
            RecipesBase.apply_recipe(KW(:sep => :red), history, :resnorm)
        )

        for data in plots
            @test length(data) == 2
            @test data[2].plotattributes[:linecolor] == :red
        end
    end
end
