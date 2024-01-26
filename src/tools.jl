# Function to split DataFrame based on column type
function _split_dataframe_by_type_(df::DataFrame)
    int_columns = [name for name in propertynames(df) if eltype(getproperty(df, name)) == Int64]
    string_columns = [name for name in propertynames(df) if eltype(getproperty(df, name)) == String]

    int_df = select(df, int_columns)
    string_df = select(df, string_columns)

    return int_df, string_df
end


# Function to accumulate number of species in one run
function _accum_loop_(data_accumulation::DataFrame)
    # Initialise count
    accum = Float64[]
    accum_species = String[]
    com = collect(1:nrow(data_accumulation))
    f(x) = sum(x) > 0

    for community_number in 1:nrow(data_accumulation)

        community = data_accumulation[community_number,:] |> DataFrame

        species  = names( # NAMES species for which abund > 0
            community[!,map(f, eachcol(community))]
        )
        append!(accum_species, species)
        accum_species = unique(accum_species)
        
        nb_accum_species = length(accum_species) # Count number of species detected so far
        
        append!(accum, nb_accum_species)
    end

    [com, accum]
end