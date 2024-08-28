mutable struct Tikz
    packages::Vector{String}
    tikzlibraries::Vector{String}
    preamble::Vector{String}
    commands::Vector{String}
    #Vertex(label::LabelType) = new(label, Arrow[], Arrow[], Dict{String,Any}())
end

#t += tikz_plot_quiver(qq,directed=true)
#add_cmd(t, tikz_plot_quiver(qq,directed=true));
#tikz_plot_quiver(t, qq, directed=true)

function tikz_plot_quiver(qq::Quiver; directed=true, vertex_color="black", linecolor="black", arrow_offset = 0.03, linewidth = 0.8, markerSize=1, draw_vertices=true ,kwargs...)
    #Arrows
    output = ""
    for arr in arrows(qq)
        if directed
            output *= tikzDrawArrow(arr.start.data["position"], arr.termination.data["position"], offset=arrow_offset, linewidth=linewidth, color="purple")
        else
            output *= tikzDrawLine( map(x->[i for i in x],collect(zip(arr.start.data["position"], arr.termination.data["position"]) ))..., color=linecolor, linewidth=linewidth)
        end
    end

    #Points
    if draw_vertices
        output *= tikzPoints(map(x -> [i for i in x],collect(zip(map(x->x.data["position"],vertices(qq))...)))..., color=vertex_color, markSize=1.6)
    end
    return output

end

function tikzDrawCommand(;linewidth=0.05, color="black", dashed=false, dashPattern=(1.0,1.0), extra="")
    output = "\\draw["
    output = output * "line width="*string(linewidth)*"pt, "
    #if dashed
    #    output = output * "dashed, "
    #end
    if dashed
        output = output * "dash pattern=on "*string(dashPattern[1])*"pt off "*string(dashPattern[2])*"pt, "
    end
    output = output * string(color) * ","
    output *= extra
    output *= "] "
    return output
end

function tikzDrawCirle(center, radius; kwargs...)
    output = tikzDrawCommand(;kwargs...) * " ("*string(center[1])*","*string(center[2])*") circle ("*string(radius)*");"
    return output
end

function tikzDrawLine(xs, ys; kwargs...)
    output = tikzDrawCommand(;kwargs...) * " "
    points = [[xs[i],ys[i]] for i in 1:length(xs)]
    for i in 1:length(points)
        output = output * "("*string(points[i][1])*","*string(points[i][2])*")"
        if i != length(points)
            output = output * " -- "
        else
            output = output * ";\n"
        end
    end
    return output
end

function tikzDrawArrow(from, to; offset=0, kwargs...)
    output = tikzDrawCommand(extra="-{Stealth[length=6pt]},"; kwargs...) * " "
    vec = (to - from);
    vec = offset*(vec/norm(vec))
    from = from + vec;
    to = to - vec
    output *= "("*string(from[1])*","*string(from[2])*")" * " -- " * "("*string(to[1])*","*string(to[2])*");\n";
    return output
end

function tikzPoint(x, y; rotate=0, color="black", marker="*", markSize=0.1, kwargs...)
    output = ""
    output = output * "\\node[mark size="*string(markSize)*"pt,color="*color*",rotate="*string(rotate)*"] at ("*string(x)*","*string(y)*") {\\pgfuseplotmark{"*marker*"}}"
    output = output * ";\n"
    return output
end

function tikzPoints(xs, ys; kwargs...)
    output = ""
    for i in 1:length(xs)
        output *= tikzPoint(xs[i], ys[i]; kwargs...)
    end
    return output
end

function compileTikzFile(tikz_filename, output_path)
    output_dir, output_filename = splitdir(output_path) # ("/../../dir", "name.tex")
    if output_dir == ""; output_dir = "."  end
    output_fullpath = joinpath(realpath(output_dir), output_filename)

    input_dir, input_filename = splitdir(tikz_filename) # ("/../../dir", "name.tex")
    if input_dir == ""; input_dir = "."  end
    input_fullpath = joinpath(realpath(input_dir), input_filename)
    input_filename_no_ext = splitext(input_filename)[1]  # ("name", "tex")

    build_dir = joinpath(tempdir(), "tikzbuild")
    cmd = `pdflatex -interaction nonstopmode -output-directory $build_dir`

    mkpath(build_dir)
    out = read(`$cmd $input_fullpath`, String)
    mv(build_dir * "/" * input_filename_no_ext * ".pdf", output_fullpath, force = true)
end

function compileTikz(tikz_string, output)
    outfile = tempname()

    open(outfile, "w") do file
        write(file, tikz_string) 
    end

    compileTikzFile(outfile, output)
end
