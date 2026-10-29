"""
    CommunityMatrix(sites, taxa, values; semantics=:abundance, metadata=Dict{String,Any}())

A validated site-by-taxon matrix with explicit measurement semantics.

Each row of `values` corresponds to one site or sampling event in `sites`,
and each column corresponds to one taxon in `taxa`. Labels are converted to
`String` and must be unique. The dimensions of `values` must match the lengths
of `sites` and `taxa`.

# Arguments
- `sites`: Site or event identifiers, in the order corresponding to the rows
  of `values`.
- `taxa`: Taxon names or identifiers, in the order corresponding to the
  columns of `values`.
- `values::AbstractMatrix{<:Real}`: Non-negative, finite measurements. Values
  are copied when the `CommunityMatrix` is constructed.
- `semantics::Symbol=:abundance`: Meaning of the matrix values. Must be one of:
  - `:incidence`: binary occurrence values; every entry must be `0` or `1`.
  - `:abundance`: abundance measurements, such as counts of individuals.
- `metadata::AbstractDict=Dict{String,Any}()`: Additional information about
  the data source or processing. Keys are converted to `String`; values are
  retained as provided.

# Throws
- `ArgumentError` if `semantics` is not supported, labels are not unique,
  any value is negative or non-finite, or `semantics` is `:incidence` and
  any value is not binary.
- `DimensionMismatch` if the matrix dimensions do not match the numbers of
  site and taxon labels.

# Notes
A zero indicates a zero measurement in the supplied matrix. Its ecological
interpretation depends on how the input data were collected; in particular,
a zero should not be treated as a verified absence unless the sampling design
supports that interpretation.
"""
struct CommunityMatrix{T<:Real,M<:AbstractMatrix{T}}
    sites::Vector{String}          # contains the identifiers of the sites or events, in the order in which they appear in the lines
    taxa::Vector{String}           # contains the names or identifiers of the taxa, in the order of the columns
    values::M                      # contains measurement values
    semantics::Symbol              # contains measurement signification (:incidence or :abundance)
    metadata::Dict{String,Any}     # contains further information on the origin or processing

    # Function to check data before creating structure
    function CommunityMatrix(
        sites, taxa, values::AbstractMatrix{T};
        semantics=:abundance,
        metadata=Dict{String,Any}()) where T<:Real
    
        # Check semantics
        semantics in (:incidence,:abundance) || throw(ArgumentError("Invalid semantics"))

        # Check dimensions: must have one row per site and one column per taxon
        size(values)==(length(sites),length(taxa)) || throw(DimensionMismatch("Labels must match matrix dimensions"))

        # Convert site labels into chains (e.g. numerical identifiers are considered as characters)
        ss,tt=String.(sites),String.(taxa)

        # Ensure that sites and taxon labels are unique
        length(unique(ss))==length(ss) && length(unique(tt))==length(tt) || throw(ArgumentError("Labels must be unique"))

        # All values must be finite and >= 0
        all(isfinite,values) && all(>=(0),values) || throw(ArgumentError("Values must be finite and non-negative"))
        
        # If incidence type data, values must be only 0 or 1
        semantics==:incidence && !all(x->x==0||x==1,values) && throw(ArgumentError("Incidence values must be binary"))
        
        new{T,typeof(copy(values))}(
            ss,tt,copy(values),semantics,
            Dict{String,Any}(string(k)=>v for (k,v) in metadata)) # keys in the metadata dictionary are converted to strings
    end
end


"""
    validate_records(table; event_col=:eventID, taxon_col=:scientificName,
                     status_col=:occurrenceStatus, count_col=nothing)

Return row-level diagnostics for a Darwin Core-style occurrence table.

The function checks that the columns used to identify sampling units and taxa
are present, and reports missing identifiers, missing or unrecognized
occurrence statuses, and invalid count values. It does not modify the input
table, remove rows, impute values, or stop execution when it finds an issue.

The input is converted to a `DataFrame`. Column arguments should be `Symbol`s
matching the input table's column names.

# Arguments
- `table`: A tabular occurrence dataset that can be converted to a
  `DataFrame`, with one row per occurrence record. The table should contain
  columns identifying the sampling unit and taxon.
- `event_col::Symbol=:eventID`: Name of the column identifying the site,
  event, or sampling unit. A missing column produces a table-level error
  diagnostic.
- `taxon_col::Symbol=:scientificName`: Name of the column identifying the
  taxon. A missing column produces a table-level error diagnostic.
- `status_col::Symbol=:occurrenceStatus`: Name of the column describing the
  occurrence status. If this column exists, non-blank values are recognized
  only when they are `"present"` or `"absent"`, ignoring case and surrounding
  whitespace. A blank status produces a warning; an unrecognized status also
  produces a warning. If the column is absent, no diagnostic is generated.
- `count_col::Union{Nothing,Symbol}=nothing`: Optional name of a column
  containing count or abundance values. If provided, the function checks that
  the column exists and that each value can be parsed as a finite,
  non-negative number. Missing, unparsable, non-finite, or negative values
  produce error diagnostics. No count values are imputed.

# Returns
A `DataFrame` with one row per diagnostic and the following columns:

- `row`: Input row number for row-level diagnostics; `missing` for diagnostics
  about absent columns.
- `column`: Name of the input column associated with the diagnostic.
- `code`: Machine-readable diagnostic code, such as `"missing_column"`,
  `"missing_identifier"`, `"missing_status"`, `"unknown_status"`, or
  `"invalid_count"`.
- `severity`: Diagnostic severity, either `"error"` or `"warning"`.
- `message`: Human-readable explanation of the diagnostic.

If no issues are detected, the function returns an empty `DataFrame` with
these columns.

# Notes
- Missing values in `event_col` or `taxon_col`, as well as blank strings, are
  reported as errors. The corresponding records are intended to be excluded
  rather than grouped under a missing identifier.
- A missing occurrence status is retained as unknown; it is not interpreted
  as an absence.
- Only the literal statuses `"present"` and `"absent"` are recognized by this
  implementation. Other Darwin Core status values are reported as
  unrecognized.
- Diagnostics are informational: this function does not itself exclude
  records. A downstream function, such as `from_records`, may apply exclusion
  rules.
- The function checks that the requested identifier and count columns exist,
  but it does not validate the full Darwin Core standard or infer suitable
  columns when they are absent.
"""
function validate_records(
    table;                          # a DarwinCore occurrence table
    event_col=:eventID,             # column name identifying the site, event or sampling unit
    taxon_col=:scientificName,      # column name identifying the taxon
    status_col=:occurrenceStatus,   # (facultative) colmn name describing the status of occurrence (present or absent)
    count_col=nothing               # (facultative) column name identifying abundance data
    )

    # Convert to DataFrame
    d=DataFrame(table);
    ns=Set(Symbol.(names(d)))

    # Create diagnostic table
    out=DataFrame(
        row=Union{Missing,Int}[],   # row number
        column=String[],            # column name
        code=String[],              # error code
        severity=String[],          # severity level
        message=String[]            # human readible explanation
    )

    # Check that the columns specified by event_col and taxon_col are present
    for c in (event_col,taxon_col)
        c in ns || push!(out,(missing,string(c),"missing_column","error","Required identifier column is absent."))
    end

    # Check that the column specified by count_col is present (if specified)
    if count_col!==nothing && !(count_col in ns)
        push!(out,(missing,string(count_col),"missing_column","error","Requested count column is absent."))
    end

    # Check rows
    for i in 1:nrow(d)

        # checks both identification columns, but only if the column exists
        for c in (event_col,taxon_col)
            c in ns && _blank(d[i,c]) && push!(out,(i,string(c),"missing_identifier","error","Missing identifiers are excluded, never grouped."))
        end
    
        # if status column included, check term consistency (must be "present" or "absent")
        if status_col in ns
            st=d[i,status_col]
            if _blank(st)
                push!(out,(i,string(status_col),"missing_status","warning","Missing status retained; absence is not inferred."))
            elseif lowercase(strip(string(st))) ∉ ("present","absent")
                push!(out,(i,string(status_col),"unknown_status","warning","Unknown status is excluded."))
            end
        end
    
        # if count column included check values
        if count_col!==nothing && count_col in ns
            x=d[i,count_col]; q=ismissing(x) ? nothing : tryparse(Float64,string(x))
            if q===nothing || !isfinite(q) || q<0
                push!(out,(i,string(count_col),"invalid_count","error","Count must be finite and non-negative; no imputation."))
            end
        end
    end
    out
end



"""
    from_records(table; event_col=:eventID, taxon_col=:scientificName,
                 status_col=:occurrenceStatus, count_col=nothing,
                 semantics=:incidence)

Build a site-by-taxon community matrix from a Darwin Core-style occurrence
table, and return the matrix together with input diagnostics and excluded
records.

The input table is converted to a `DataFrame`. Each retained row contributes
to the value for its event–taxon pair. The resulting `CommunityMatrix` records
the selected measurement semantics in its `semantics` field.

Zeros in the resulting matrix mean that no retained occurrence record was
found for the corresponding event–taxon pair. They do **not** represent
verified absences.

# Arguments
- `table`: A tabular occurrence dataset convertible to a `DataFrame`, usually
  with one row per occurrence record. It must contain the columns specified
  by `event_col` and `taxon_col`.
- `event_col::Symbol=:eventID`: Name of the column identifying the site,
  event, or sampling unit. Its non-blank values become the rows of the
  community matrix.
- `taxon_col::Symbol=:scientificName`: Name of the column identifying the
  taxon. Its non-blank values become the columns of the community matrix.
- `status_col::Symbol=:occurrenceStatus`: Optional name of the column
  describing occurrence status. If the column exists, `"present"` records
  are eligible for inclusion and `"absent"` records are excluded. Status
  matching ignores letter case and surrounding whitespace. Blank or missing
  statuses are retained; they are not interpreted as absences. Other
  non-blank statuses are excluded as unknown. If this column is absent, rows
  with valid event and taxon identifiers are processed without a status
  filter.
- `count_col::Union{Nothing,Symbol}=nothing`: Optional name of the abundance
  column. It is required when `semantics=:abundance`. In that mode, each
  retained value must parse as a finite, non-negative number; values for the
  same event–taxon pair are summed.
- `semantics::Symbol=:incidence`: Measurement represented in the community
  matrix. Must be either `:incidence` or `:abundance`.
  - `:incidence`: Each event–taxon pair with one or more retained records is
    assigned `1.0`, regardless of the number of records.
  - `:abundance`: The values in `count_col` are summed for each event–taxon
    pair.

# Returns
A 3-tuple:

1. `community::CommunityMatrix`: The site-by-taxon matrix. Site and taxon
   labels are converted to strings and have leading and trailing whitespace
   removed. Labels are ordered by their first appearance among retained
   records. The matrix values are floating-point numbers.
2. `diagnostics::DataFrame`: Row-level and table-level diagnostics returned
   by [`validate_records`](@ref), based on the input table and the specified
   column arguments. Diagnostics report issues; they do not themselves
   determine which records are excluded.
3. `exclusions::DataFrame`: Records omitted while constructing the matrix,
   with columns `source_row` (the row number in the input table) and `reason`
   (the exclusion reason). Records can be excluded for missing identifiers,
   explicit absence, unknown status, or invalid or missing abundance values.

# Behavior
- A missing `event_col` or `taxon_col` raises an `ArgumentError`.
- If `semantics` is not `:incidence` or `:abundance`, an `ArgumentError` is
  raised.
- `semantics=:abundance` requires `count_col` to be specified and present in
  the table; otherwise, an `ArgumentError` is raised.
- Rows with a blank or missing event or taxon identifier are excluded and
  are not grouped under a shared missing-value label.
- Rows with status `"absent"` or an unrecognized non-blank status are
  excluded. A missing or blank status is retained, because absence must not
  be inferred from a missing status.
- In abundance mode, a row with a missing, unparsable, non-finite, or
  negative count is excluded. No values are imputed.
- Multiple retained records for the same event–taxon pair are collapsed into
  one matrix cell. Incidence mode records whether any such record exists;
  abundance mode sums their counts.
- The returned matrix metadata includes the source format, the interpretation
  of zero values, the number of input rows, and the number of excluded rows.

# Examples
```julia
community, diagnostics, exclusions = from_records(occurrences)

community, diagnostics, exclusions = from_records(
    occurrences;
    event_col=:eventID,
    taxon_col=:scientificName,
    semantics=:abundance,
    count_col=:individualCount,
)
"""
function from_records(
    table;                          # a DarwinCore occurrence table
    event_col=:eventID,             # column name identifying the site, event or sampling unit
    taxon_col=:scientificName,      # column name identifying the taxon
    status_col=:occurrenceStatus,   # (facultative) colmn name describing the status of occurrence (present or absent)
    count_col=nothing,              # (facultative) column name identifying abundance data
    semantics=:incidence            # measurement signification (:incidence or :abundance)
    )

    # Check semantics and convert to table
    semantics in (:incidence,:abundance)||throw(ArgumentError("Invalid semantics"))
    d=DataFrame(table); ns=Set(Symbol.(names(d)))

    # Apply diagnostic function
    diag=validate_records(d;event_col,taxon_col,status_col,count_col)

    # check presence of required columns
    all(c->c in ns,(event_col,taxon_col))||throw(ArgumentError("Required identifier columns absent"))

    # if abundance analysis are chosen, check required column
    semantics==:abundance && (count_col===nothing || !(count_col in ns)) && throw(ArgumentError("Abundance requires an explicit count column"))

    # Convert to useful types
    sites=String[];
    taxa=String[]; 
    vals=Dict{Tuple{String,String},Float64}()
    exc=DataFrame(source_row=Int[],reason=String[])

    """
    Goes through each line and adds the line to the exclusions if:

    - its event ID or taxon ID is empty;
    - its status is ‘absent’;
    - its status is neither ‘present’ nor ‘absent’;
    - in ‘abundance’ mode, its count is missing, non-numeric, infinite or negative.
    """
    for i in 1:nrow(d)
        ev,tx=d[i,event_col],d[i,taxon_col]
        if _blank(ev)||_blank(tx); push!(exc,(i,"missing identifier"));continue;end
        st=status_col in ns ? d[i,status_col] : missing
        if !_blank(st)
            status=lowercase(strip(string(st)))
            if status=="absent";push!(exc,(i,"explicit absence"));continue;end
            if status!="present";push!(exc,(i,"unknown status"));continue;end
        end
        v=1.0
        if semantics==:abundance
            raw=d[i,count_col]; q=ismissing(raw) ? nothing : tryparse(Float64,string(raw))
            if q===nothing||!isfinite(q)||q<0;push!(exc,(i,"invalid or missing count"));continue;end
            v=q
        end
        e,t=strip(string(ev)),strip(string(tx))
        e in sites||push!(sites,e);t in taxa||push!(taxa,t)
        key=(e,t);vals[key]=get(vals,key,0.0)+v
    end

    # Aggregate cleaned records
    m=zeros(length(sites),length(taxa));si=Dict(x=>i for(i,x) in enumerate(sites));ti=Dict(x=>i for(i,x) in enumerate(taxa))

    for ((e,t),v) in vals;m[si[e],ti[t]]=semantics==:incidence ? 1.0 : v;end

    # Build metadata matrix
    meta=Dict(
        "source_format"=>"Darwin Core-style table",
        "zero_interpretation"=>"unrecorded, not verified absence",
        "input_rows"=>nrow(d),
        "excluded_rows"=>nrow(exc)
    )

    CommunityMatrix(sites, taxa, m; semantics,metadata=meta),
    diag,
    exc

end

"""
    from_table(path; delim='\\t', kwargs...)

Read a tab-delimited text file and build a community matrix from its
occurrence records.

The file is loaded into a `DataFrame` and passed to [`from_records`](@ref).
This function supports tab-separated text files; it does not read complete
Darwin Core Archive (DwC-A) archives.

# Arguments
- `path`: Path to the input `.txt` file. The file must exist.
- `delim='\\t'`: Field delimiter passed to `CSV.read`. By default, fields are
  separated by tab characters.
- `kwargs...`: Keyword arguments forwarded to [`from_records`](@ref), such as
  `event_col`, `taxon_col`, `status_col`, `count_col`, and `semantics`.

# Returns
Returns the same 3-tuple as [`from_records`](@ref):
1. A `CommunityMatrix`.
2. A diagnostics `DataFrame`.
3. An exclusions `DataFrame`.

# Errors
- Throws an `ArgumentError` if `path` does not refer to an existing file.
- May throw errors from `CSV.read` if the file cannot be read or parsed.
- May throw errors from `from_records` if its required columns or arguments
  are invalid.

# Example
```julia
community, diagnostics, exclusions = from_table(
    "occurrences.txt";
    event_col=:eventID,
    taxon_col=:scientificName,
    semantics=:abundance,
    count_col=:individualCount,
)
"""

function from_table(path; delim='\t', kwargs...)
    isfile(path) || throw(ArgumentError("File not found"))
    from_records(CSV.read(path, DataFrame; delim=delim); kwargs...)
end