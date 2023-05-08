using ChEMBL

# Load the Tox21 dataset
tox21_data = CSV.File("tox21.csv", header=true) |> DataFrame

# Extract the ChEMBL IDs for each compound
chembl_ids = []
for smiles in tox21_data[:smiles]
    result = search_compounds(smiles, search_type="similarity", threshold=0.9)
    if length(result) > 0
        push!(chembl_ids, result[1]["chembl_id"])
    else
        push!(chembl_ids, "")
    end
end
