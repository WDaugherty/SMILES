using Flux
using Statistics
using DataFrames
using Random

# Load the preprocessed data
tox21_data = CSV.File("tox21_data/tox21_10k_data_with_descriptors.csv", header=true) |> DataFrame

# Define the input and output variables for the neural network
inputs = Matrix(tox21_data[[:$(name) for name in keys(descriptors)]] |> Matrix)
outputs = Matrix(tox21_data[:, [:$(name) for name in ["NR.AhR", "NR.AR", "NR.AR.LBD", "NR.Aromatase", "NR.ER", "NR.ER.LBD", "NR.PPAR.gamma", "SR.ARE", "SR.ATAD5", "SR.HSE", "SR.MMP", "SR.p53"]]])

# Define the neural network architecture
n_inputs = size(inputs, 2)
n_hidden = 50
n_outputs = size(outputs, 2)
model = Chain(
    Dense(n_inputs, n_hidden, relu),
    Dropout(0.5),
    Dense(n_hidden, n_hidden, relu),
    Dropout(0.5),
    Dense(n_hidden, n_outputs, sigmoid)
)

# Define the loss function
function binary_cross_entropy(y_pred, y_true)
    -mean(y_true .* log.(y_pred) .+ (1 .- y_true) .* log.(1 .- y_pred))
end

# Define the optimizer and the training function
optimizer = ADAM()
function train_model(model, inputs, outputs, n_epochs=100, batch_size=128)
    for epoch in 1:n_epochs
        # Shuffle the data
        indices = shuffle(1:size(inputs, 1))
        # Train the model in mini-batches
        for batch_start in 1:batch_size:size(inputs, 1)-batch_size+1
            batch_end = batch_start+batch_size-1
            batch_indices = indices[batch_start:batch_end]
            x_batch = inputs[batch_indices, :]
            y_batch = outputs[batch_indices, :]
            loss = binary_cross_entropy(model(x_batch), y_batch)
            Flux.back!(loss)
            Flux.update!(optimizer, model)
        end
        # Print the loss after each epoch
        println("Epoch $epoch: Loss $(binary_cross_entropy(model(inputs), outputs))")
    end
end

# Train the neural network
train_model(model, inputs, outputs, n_epochs=100)
