function CountSpheroids(data::DataFrame)

    CellLines         = unique(data.CellLine)
    InitialConditions = sort(unique(data.InitialCondition))
    Days              = sort(unique(data.Day))

    df = DataFrame(
        CellLine         = vcat([fill(CellLine,length(InitialConditions)) for CellLine ∈ CellLines]...),
        InitialCondition = [InitialConditions; InitialConditions]
    )

    for day ∈ Days
        df[:,Symbol("D$day")] = [nrow(@subset(data, :CellLine .== df.CellLine[i], :InitialCondition .== df.InitialCondition[i], :Day .== day)) for i = 1:nrow(df)]
    end

    df

end