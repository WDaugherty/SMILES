using CSV
using Plots

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Load the chemical data
molweights = [430.5, 356.9, 405.5, 505.6, 315.8, 401.4, 494.5, 412.6, 383.5, 305.4]
boilingpoints = [203.5, 259.2, 237.5, 402.9, 63.5, 364.6, 601.4, 387.4, 319.8, 210.9]

# Create the heatmap
heatmap(["MolWeight", "BoilingPoint"], tox21_data[:id], [molweights boilingpoints], color=:viridis, xlabel="", ylabel="Compound ID", title="Correlation between chemical features and toxicity levels")
