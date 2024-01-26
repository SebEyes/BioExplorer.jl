# Function to split DataFrame based on column type
function split_dataframe_by_type(df::DataFrame)
    int_columns = [name for name in propertynames(df) if eltype(getproperty(df, name)) == Int64]
    string_columns = [name for name in propertynames(df) if eltype(getproperty(df, name)) == String]

    int_df = select(df, int_columns)
    string_df = select(df, string_columns)

    return int_df, string_df
end