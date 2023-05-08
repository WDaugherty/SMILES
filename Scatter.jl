using CSV
using Plots

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Load the chemical data
molweights = [430.5, 356.9, 405.5, 505.6, 315.8, 401.4, 494.5, 412.6, 383.5, 305.4]

# Create the scatter plot
scatter(molweights, tox21_data[:toxicity], xlabel="Molecular weight", ylabel="Toxicity level")
