#=
    Table 1

    Summarise number of spheroids (both all and sub-sampled)

=#

using Latexify

confocal = CSV.read("Data/ConfocalData.csv",DataFrame)
confocal_all = CSV.read("Data/AllConfocalData.csv",DataFrame)

include("../../Data/CountSpheroids.jl")

counts = CountSpheroids(confocal)
counts_all = CountSpheroids(confocal_all)

tabS1 = [counts_all; counts]
tabS1.Totals = sum(Array(tabS1[:,3:end]),dims=2)[:]

write("$(@__DIR__)/tabS1.tex",latexify(tabS1,env=:tabular))