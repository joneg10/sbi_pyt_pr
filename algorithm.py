import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import torch.optim as optim
from sklearn.model_selection import train_test_split
from read_pdb_classes import *

data = pd.read_csv('../output.csv' , na_values = 0 , sep = "," , )
data.fillna(0, inplace=True)

# From TrainingSet class to format the training set from pdb files, still to computational expensive to run in the server:
#paths = TrainingSet("./pdb_ids")
#training_set = paths.get_formated_set()

# Set "is_lbs" as the last column
data["is_lbs"] = data.pop("is_lbs")

# Set "Unnamed: 0" as row names
data = data.set_index("Unnamed: 0")

# Scale the data
data.loc[:, data.columns != "is_lbs"] = (data.loc[:, data.columns != "is_lbs"] - data.loc[:, data.columns != "is_lbs"].mean()) / data.loc[:, data.columns != "is_lbs"].std()

### Neural network ###

# Define the model

num_data = np.array(data)
X = num_data[:, :-1]
Y = num_data[:, -1]

X = torch.tensor(X, dtype=torch.float32)
Y = torch.tensor(Y, dtype=torch.float32).reshape(-1, 1)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, random_state=42)


model = nn.Sequential(
    nn.Linear(46, 69),
    nn.ReLU(),
    nn.Linear(69, 46),
    nn.ReLU(),
    nn.Linear(46, 1),
    nn.Sigmoid()
)


loss_fn = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)

# training

n_epochs = 100
batch_size = 10
 
for epoch in range(n_epochs):
    for i in range(0, len(X_train), batch_size):
        Xbatch = X_train[i:i+batch_size]
        Y_pred = model(Xbatch)
        Ybatch = Y_train[i:i+batch_size]
        loss = loss_fn(Y_pred, Ybatch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    print(f'Finished epoch {epoch}, latest loss {loss}')


with torch.no_grad():
    Y_pred = model(X_test)
 
accuracy = (Y_pred.round() == Y_test).float().mean()
print(f"Accuracy {accuracy}")

model.state_dict()

torch.save(model.state_dict(), "neural_network.pytorch")