cd(@__DIR__)

using Pkg
Pkg.activate(".")

using Printf

##################################
# Data Structures & Constructors #
##################################

struct GroLine
    resNr::Int
    resNm::String
    atmNm::String
    atmNr::Int
    pos::Tuple{Float32,Float32,Float32}
end

function GroLine(line::String)
    resNr = parse(Int, line[1:5])
    resNm = strip(line[6:10])
    atmNm = strip(line[11:15])
    atmNr = parse(Int, line[16:20])
    x = parse(Float32, line[21:28])
    y = parse(Float32, line[29:36])
    z = parse(Float32, line[37:44])
    return GroLine(resNr,resNm,atmNm,atmNr,(x,y,z))
end

struct Resid
    resNm::String
    resNr::Int
    pos::Tuple{Float32, Float32, Float32}
    mass::Float32
end

#############################
# Functions to Read & Write #
#############################

function readGroFile(fileName::String)
    groLines = GroLine[]
    open(fileName, "r") do input
        readline(input) # skip the title
        nrAtms = parse(Int, readline(input))
        for _ in 1:nrAtms
            push!(groLines, GroLine(readline(input)))
        end
    end
    return groLines
end

function writeMasses(out::IOStream,atmTypes::Vector{Tuple{String,Float32}})
    println(out, "Masses\n")
    for (i,t) in enumerate(atmTypes)
        @printf(out, "%6d %10.3f   #%s\n", i, t[2], t[1])
    end
    println(out, "")
end

function writeSimulBox(out::IOStream,simulBox)
    @printf(out, "%8.1f %8.1f %s\n", simulBox.xlo, simulBox.xhi, "  xlo  xhi")
    @printf(out, "%8.1f %8.1f %s\n", simulBox.ylo, simulBox.yhi, "  ylo  yhi")
    @printf(out, "%8.1f %8.1f %s\n", simulBox.zlo, simulBox.zhi, "  zlo  zhi") 
    println(out, "")   
end

###################################
# Functions to build the CG model #
###################################

function findCA(fileName::String)
    masses = Dict([("H",1.0080f0),("C",12.011f0),("N",14.007f0),("O",15.999f0),("S",32.06f0)])
    out = Resid[] 
    groLines = readGroFile(fileName)
    res0 = -1; m0 = 0.0f0; nm0 = ""; pos0 = nothing
    for groLine in groLines
        if groLine.resNr != res0
            isnothing(pos0) || push!(out, Resid(nm0,res0,pos0,m0))
            m0 = groLine.atmNm[1:1] in keys(masses) ? masses[groLine.atmNm[1:1]] : 0.0f0
            res0 = groLine.resNr; nm0 = groLine.resNm; pos0 = nothing
        else
            m0 += groLine.atmNm[1:1] in keys(masses) ? masses[groLine.atmNm[1:1]] : 0.0f0
        end
        (groLine.atmNm == "CA") && (pos0 = groLine.pos)
    end
    isnothing(pos0) || push!(out, Resid(nm0,res0,pos0,m0))
    return out
end

function fixAtomTypes(residueList::Vector{Resid})
    atmTypes = Set{Tuple{String,Float32}}()
    for res in residueList
        push!(atmTypes, (res.resNm, res.mass))
    end
    return sort(collect(atmTypes))
end

function fixSimulBox(residueList::Vector{Resid}; cutoff=10.0)
    xv = [res.pos[1] for res in residueList]
    yv = [res.pos[2] for res in residueList]
    zv = [res.pos[3] for res in residueList]
    # NB: positions in .gro are in nm, but LAMMPS uses Å
    xmin = minimum(xv)*10; xmax = maximum(xv)*10; xave = (xmin+xmax)/2 
    ymin = minimum(yv)*10; ymax = maximum(yv)*10; yave = (ymin+ymax)/2 
    zmin = minimum(zv)*10; zmax = maximum(zv)*10; zave = (zmin+zmax)/2
    l = sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2)+2*cutoff
    xlo = xave-l/2; ylo = yave-l/2; zlo = zave-l/2
    xhi = xave+l/2; yhi = yave+l/2; zhi = zave+l/2
    return (xlo=xlo,xhi=xhi,ylo=ylo,yhi=yhi,zlo=zlo,zhi=zhi)
end

function createCGmodel(frm0, outputFile)
    resLst = findCA(frm0)
    simulBox = fixSimulBox(resLst)
    atmTypes = fixAtomTypes(resLst)
    open(outputFile, "w") do output
        writeSimulBox(output, simulBox)
        writeMasses(output, atmTypes)
    end
end