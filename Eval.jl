using ROCAnalysis
using DataFrames
using Plots

# Load the preprocessed data
tox21_data = CSV.File("tox21_data/tox21_10k_data_with_descriptors.csv", header=true) |> DataFrame

# Define the input and output variables for the neural network
inputs = Matrix(tox21_data[[:$(name) for name in keys(descriptors)]] |> Matrix)
outputs = Matrix(tox21_data[:, [:$(name) for name in ["NR.AhR", "NR.AR", "NR.AR.LBD", "NR.Aromatase", "NR.ER", "NR.ER.LBD", "NR.PPAR.gamma", "SR.ARE", "SR.ATAD5", "SR.HSE", "SR.MMP", "SR.p53"]]])

# Evaluate the model on the test data
y_pred = model(inputs)
auc_scores = [rocauc(outputs[:, i], y_pred[:, i]) for i in 1:size(outputs, 2)]

# Print the AUC scores for each toxicity endpoint
for (i, endpoint) in enumerate(["NR.AhR", "NR.AR", "NR.AR.LBD", "NR.Aromatase", "NR.ER", "NR.ER.LBD", "NR.PPAR.gamma", "SR.ARE", "SR.ATAD5", "SR.HSE", "SR.MMP", "SR.p53"])
    println("AUC for $endpoint: $(auc_scores[i])")
end

# Plot the ROC curves for each toxicity endpoint
roc_curves = [ROCAnalysis.roc(outputs[:, i], y_pred[:, i]) for i in 1:size(outputs, 2)]
plot(roc_curves, layout=(4, 3), size=(800, 600))
