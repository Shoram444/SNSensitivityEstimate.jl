using SNSensitivityEstimate
using Test

@testset "SNSensitivityEstimate.jl" begin
    # Write your tests here.
    SNparams_test = SNSensitivityEstimate.SNparams
    @test SNparams_test["Nₐ"] == 6.02214e23
end
