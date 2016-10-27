using Base.Test

# notebooks
@testset "example notebooks" begin
    using IJulia

    outputdir = joinpath(tempdir(), "VariableHeightInvertedPendulum")
    if !isdir(outputdir)
        mkpath(outputdir)
    end
    jupyter = IJulia.jupyter
    notebookdir = joinpath("..", "notebook")
    for f in filter(x -> endswith(x, "ipynb"), readdir(notebookdir))
        notebook = joinpath(notebookdir, f)
        output = joinpath(outputdir, f)
        # skip on nightly because notebooks specify version 0.5
        @test begin run(`$jupyter nbconvert --to notebook --execute $notebook --output $output --ExecutePreprocessor.timeout=90`); true end
    end
end
