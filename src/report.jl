"""
    write_report(r::AnalysisResult, path;
                 title="BioExplorer community analysis",
                 overwrite=false)

Render an `AnalysisResult` as a PDF report using Typstry.jl and Typst.

Each table is limited to its first 35 rows and first 8 columns to preserve
the report layout. A notice is added when rows or columns are omitted.
The tables in `r` are not modified.

# Arguments
- `r`: Analysis results to include in the report.
- `path`: Destination PDF path; the filename must end in `.pdf`.
- `title`: Title displayed at the beginning of the report.
- `overwrite`: If `false`, refuse to write if the destination already exists.

# Returns
The destination path.

# Throws
- `ArgumentError` if `path` does not end in `.pdf`.
- `ArgumentError` if the destination exists and `overwrite` is `false`.
- An error if Typst compilation fails.
"""
function write_report(
    r::AnalysisResult,
    path;
    title = "BioExplorer community analysis",
    overwrite = false,
)
    endswith(lowercase(path), ".pdf") ||
        throw(ArgumentError("Destination must be a PDF"))

    ispath(path) && !overwrite &&
        throw(ArgumentError("Refusing to overwrite destination"))

    # Escape text before inserting it into Typst source.
    typst_literal(x) =
        "\"" * replace(
            string(x),
            "\\" => "\\\\",
            "\"" => "\\\"",
            "\n" => "\\n",
            "\r" => "\\r",
            "\t" => "\\t",
        ) * "\""

    function write_section(io, section_title, table)
        println(io, "== ", typst_literal(section_title))

        if ncol(table) == 0 || nrow(table) == 0
            println(io, "No results available.\n")
            return
        end

        max_rows = 35
        max_columns = 8
        displayed_rows = min(nrow(table), max_rows)
        displayed_columns = min(ncol(table), max_columns)

        # Build a table using only the selected rows and columns.
        println(io, "#table(columns: ", displayed_columns, ", inset: 4pt,")

        for column_index in 1:displayed_columns
            println(io, "text(", typst_literal(names(table)[column_index]), "),")
        end

        for row_index in 1:displayed_rows
            for column_index in 1:displayed_columns
                value = table[row_index, column_index]
                display_value =
                    ismissing(value) ? "NA" :
                    value isa AbstractFloat ? string(round(value; digits = 4)) :
                    string(value)

                println(io, "text(", typst_literal(display_value), "),")
            end
        end

        println(io, ")\n")

        if nrow(table) > max_rows || ncol(table) > max_columns
            println(
                io,
                "*Table truncated: showing ",
                displayed_rows,
                " of ",
                nrow(table),
                " rows and ",
                displayed_columns,
                " of ",
                ncol(table),
                " columns.*\n",
            )
        end
    end

    # Assemble the Typst document.
    source = sprint() do io
        println(io, "#set page(paper: \"a4\", flipped: true, margin: 1.5cm)")
        println(io, "#set text(size: 9pt)")
        println(io, "= #text(", typst_literal(title), ")")
        println(
            io,
            "IN DEVELOPMENT: independently validate methods and input sampling design before interpretation.",
        )
        println(
            io,
            "Zeros represent unrecorded taxa, not verified absence. ",
            "Accumulation SD is not a confidence interval.",
        )

        write_section(io, "Alpha diversity", r.alpha)
        write_section(io, "Beta dissimilarity", r.beta)
        write_section(io, "Taxon composition", r.composition)
        write_section(io, "Site accumulation", r.accumulation)

        if r.estimators !== nothing
            write_section(io, "Richness estimators", r.estimators)
        end

        write_section(io, "Input diagnostics", r.diagnostics)
    end

    destination = abspath(path)
    mkpath(dirname(destination))

    # Typstry's @run invokes Typst_jll-backed TypstCommand.
    mktempdir() do temporary_directory
        source_path = joinpath(temporary_directory, "report.typ")
        pdf_path = joinpath(temporary_directory, "report.pdf")

        write(source_path, source)

        run(Typstry.TypstCommand(["compile", source_path, pdf_path]))
        
        isfile(pdf_path) ||
            error("Typst did not produce the expected PDF: $pdf_path")

        cp(pdf_path, destination; force = overwrite)
    end

    return path
end