using RDKit
using DataFrames
using CSV
using Flux
using Statistics
using Plots
using HTTP

# Load the pre-trained model
pretrained_model = Flux.load("pretrained_model.jld2")["model"]

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Convert the SMILES strings to molecular graphs
molgraphs = [smiles_to_molgraph(smiles) for smiles in tox21_data[:smiles]]

# Extract the molecular descriptors
descriptors = [extract_descriptors(molgraph) for molgraph in molgraphs]

# Convert the descriptors to a matrix
inputs = Matrix([values(desc) for desc in descriptors])

# Compute the feature importance scores using integrated gradients
function integrated_gradients(model, inputs, outputs, baseline=zeros(size(inputs, 2)), n_steps=50)
    gradient_fn = gradient(() -> model(inputs), inputs)[1]
    interpolated_inputs = [baseline .+ (i/n_steps)*(inputs .- baseline) for i in 0:n_steps]
    interpolated_outputs = model.(interpolated_inputs)
    gradient_outputs = gradient_fn.(interpolated_inputs)
    return (inputs .- baseline) .* mean(gradient_outputs, dims=1)
end

# Compute the feature importance scores for the Tox21 dataset
scores = zeros(size(inputs))
for i in 1:size(outputs, 2)
    scores[:, i] = integrated_gradients(pretrained_model, inputs, outputs[:, i], baseline=zeros(size(inputs, 2)))
end

# Plot the feature importance scores
heatmap([keys(descriptors)'; scores'], aspect_ratio=1, xaxis=:auto, yaxis=:auto, color=:viridis, xlabel="Molecular Descriptors", ylabel="Toxicity Endpoints", title="Integrated Gradients Feature Importance Scores")
