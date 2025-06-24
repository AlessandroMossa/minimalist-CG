cd(@__DIR__)

using Pkg
Pkg.activate(".")

using Printf

##################################
# Data Structures & Constructors #
##################################

struct Conf
    reference_structure::String
    box_margin::Float32
    trans_peptide_bond::Vector{Float32}
    cis_peptide_bond::Vector{Float32}
end

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
    # a .gro file is written in nm, but LAMMPS works in Å
    x = parse(Float32, line[21:28])*10
    y = parse(Float32, line[29:36])*10
    z = parse(Float32, line[37:44])*10
    return GroLine(resNr,resNm,atmNm,atmNr,(x,y,z))
end

struct Resid
    resNm::String
    resNr::Int
    pos::Tuple{Float32, Float32, Float32}
    mass::Float32
end

####################
# Vector Functions #
####################

⋅(u,v) = u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
×(u,v) = (u[2]*v[3]-v[2]*u[3], u[3]*v[1]-v[3]*u[1], u[1]*v[2]-v[1]*u[2])
norm(v) = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
angle(u,v) = acosd(clamp(u⋅v/(norm(u)*norm(v)), -1, 1))

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

function writeBondCoeffs(out::IOStream,conf::Conf)
    println(out, "Bond Coeffs\n")
    @printf(out, "     1 harmonic %8.3f %8.3f\n", conf.trans_peptide_bond...)
    @printf(out, "     2 harmonic %8.3f %8.3f\n", conf.cis_peptide_bond...)
    println(out, "")
end

function writeAngleCoeffs(out::IOStream,θrange::AbstractRange{Float32})
    println(out, "Angle Coeffs    # cosine/squared\n")
    i = 0
    for θ in θrange
        i += 1
        @printf(out, "%5d %8.2f %8.1f\n", i, kθ(θ/180*π), θ)
    end
    println(out, "")
end

function writeAtoms(out::IOStream,residueList::Vector{Resid},atmTypes::Vector{Tuple{String,Float32}})
    println(out, "Atoms\n")
    for (i,res) in enumerate(residueList)
        j = only(indexin([(res.resNm,res.mass)], atmTypes))
        @printf(out, "%5d  1 %5d %8.3f %8.3f %8.3f  0  0  0\n", i, j, res.pos...)
    end
    println(out,"")
end

function writeBonds(out::IOStream,residueList::Vector{Resid},conf::Conf)
    thr = (conf.trans_peptide_bond[2] + conf.cis_peptide_bond[2])/2
    println(out, "Bonds\n")
    for i in 1:length(residueList)-1
        l = norm(residueList[i+1].pos .- residueList[i].pos)
        j = (l > thr) ? 1 : 2
        @printf(out, "%5d %4d %5d %5d\n", i, j, i, i+1)
    end
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

function fixSimulBox(residueList::Vector{Resid}; box_margin=10.0)
    xv = [res.pos[1] for res in residueList]
    yv = [res.pos[2] for res in residueList]
    zv = [res.pos[3] for res in residueList]
    xmin = minimum(xv); xmax = maximum(xv); xave = (xmin+xmax)/2 
    ymin = minimum(yv); ymax = maximum(yv); yave = (ymin+ymax)/2 
    zmin = minimum(zv); zmax = maximum(zv); zave = (zmin+zmax)/2
    l = sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2)+2*box_margin
    xlo = xave-l/2; ylo = yave-l/2; zlo = zave-l/2
    xhi = xave+l/2; yhi = yave+l/2; zhi = zave+l/2
    return (xlo=xlo,xhi=xhi,ylo=ylo,yhi=yhi,zlo=zlo,zhi=zhi)
end

function fixThetaAngles(residueList::Vector{Resid})
    θv = Float32[]
    for i in 1:length(residueList)-2
        push!(θv, angle(residueList[i].pos .- residueList[i+1].pos, 
                        residueList[i+2].pos .- residueList[i+1].pos ))
    end
    return θv
end

function kθ(θ)
    B = 3000 # kcal/mol
    k0 = 10  # kcal/mol
    β = 5/3  
    return (k0 + B * (sin(β*θ)/(β*θ))^2)/(2*sin(θ)^2)
end

function prepareConfiguration()
    reference_structure = "md_0.gro"
    box_margin =          10.0f0          #Å
    trans_peptide_bond = [98.5f0, 3.80f0] #Å, kcal/mol
    cis_peptide_bond =   [70.0f0, 2.88f0] #Å, kcal/mol
    return Conf(reference_structure,box_margin,trans_peptide_bond,cis_peptide_bond)
end

function createCGmodel(outputFile::String)
    conf = prepareConfiguration()
    resLst = findCA(conf.reference_structure)
    simulBox = fixSimulBox(resLst; box_margin=conf.box_margin)
    atmTypes = fixAtomTypes(resLst)
    θAngles = fixThetaAngles(resLst)
    θrange = (floor(minimum(θAngles))+0.5f0):1.0f0:(ceil(maximum(θAngles))+0.5f0)
    open(outputFile, "w") do output
        writeSimulBox(output, simulBox)
        writeMasses(output, atmTypes)
        writeBondCoeffs(output, conf)
        writeAngleCoeffs(output, θrange)
        writeAtoms(output,resLst,atmTypes)
        writeBonds(output,resLst,conf)

    end
end
