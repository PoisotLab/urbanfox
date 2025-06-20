using SpeciesDistributionToolkit
using BiodiversityObservationNetworks
using CairoMakie
using Statistics
const SDT = SpeciesDistributionToolkit

# Get some polygon data
berlin = getpolygon(PolygonData(OpenStreetMap, Places), place="Berlin")
bbox = SDT.boundingbox(berlin, padding=0.05)

# Read the pharos data
if ~isfile("data/occurrences.json")
    include("helpers/cleandata.jl")
end
occ = OccurrencesInterface.load("data/occurrences.json")

# Reality check
fig = Figure()
ax = Axis(fig[1,1]; aspect=1.3)
lines!(ax, berlin, color=:black)
scatter!(ax, occ, color=presence(occ))
fig

# Grid for k-means
prevalence = SDMLayer(zeros(Float64, 100, 100), x=(bbox.left, bbox.right), y=(bbox.bottom, bbox.top))
mask!(prevalence, berlin)

for k in keys(prevalence)
    x, y = eastings(prevalence)[k[2]], northings(prevalence)[k[1]]
    # Pick closest neighbors
    d = [PseudoAbsences._distancefunction((x, y), p) for p in place(occ)]
    # Get five closest neighbors
    closest = partialsortperm(d, 1:40)
    # Outcomes, distances
    out = occ[closest]
    dst = 1.0 ./ d[closest]
    # Get prevalence weighted by inverse distance
    prevalence[k] = sum((presence(occ[closest]) .* d[closest]) ./ sum(d[closest]))
end

# z-score for prevalence
z = (prevalence .- mean(prevalence)) ./ std(prevalence)

fig = Figure()
ax = Axis(fig[1,1])
hm = heatmap!(ax, z, colormap=:diverging_gwv_55_95_c39_n256, colorrange=(-2, 2))
lines!(ax, berlin, color=:black)
hidespines!(ax)
hidedecorations!(ax)
Colorbar(fig[1,2], hm)
fig

# Do a BON
num_nodes = 50
bon = sample(CubeSampling(num_nodes), [Raster(z)])