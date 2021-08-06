#=
    SubsetData.jl
=#

using Distributions, StatsBase, Random
include("CountSpheroids.jl")

confocal_all = CSV.read("Data/AllConfocalData.csv",DataFrame)


# Maximum 10 each day up to day 18
maxdata = (CellLine,Day,InitialCondition) -> ((CellLine == "983b") & (Day == 18) & (InitialCondition ≥ 5000))  | 
                                             ((CellLine == "983b") & (Day == 21) & (InitialCondition == 2500)) | 
                                             ((CellLine == "793b") & (Day == 24)) ? 20 : 10

# Use reproducable random numbers
rng = MersenneTwister(0)

# Randomly remove rows
confocal = DataFrame()
for CellLine ∈ unique(confocal_all.CellLine), Day ∈ sort(unique(confocal_all.Day)), InitialCondition ∈ sort(unique(confocal_all.InitialCondition))

    global confocal

    confocal_sub = @subset(confocal_all, :CellLine .== CellLine, :Day .== Day, :InitialCondition .== InitialCondition)
    nrows = nrow(confocal_sub)

    if nrows > maxdata(CellLine,Day,InitialCondition)
        confocal_sub = confocal_sub[sample(rng,1:nrows,maxdata(CellLine,Day,InitialCondition),replace=false),:]
    end

    confocal = vcat(confocal,confocal_sub)

end

CSV.write("Data/ConfocalData.csv",confocal)

counts_full   = CountSpheroids(confocal_all)
counts_subset = CountSpheroids(confocal)
