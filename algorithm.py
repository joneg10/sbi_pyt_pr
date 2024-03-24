import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import torch.optim as optim
from sklearn.model_selection import train_test_split
from create_training_set import TrainingSet
from sklearn.preprocessing import StandardScaler
from torch.optim.lr_scheduler import ExponentialLR
from sklearn.metrics import roc_curve, roc_auc_score, auc
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt




# data = pd.read_csv('../output.csv' , na_values = 0 , sep = "," , )
# data.fillna(0, inplace=True)

# From TrainingSet class to format the training set from pdb files, still too computational expensive to run in the server:
paths = TrainingSet("../scpdb_files_1000")

training_set = paths.get_formated_set()
training_set.fillna(0, inplace=True)

# I wanted to save it to parquet to avoid running the previous code again
output_file = '../trainingSet_20240324.parquet'
training_set.to_parquet(output_file)

# Set "is_lbs" as the last column
training_set["is_lbs"] = training_set.pop("is_lbs")

training_set = pd.read_parquet('./trainingSet.parquet')
### Neural network ###

# Define the model

num_training_set = np.array(training_set)
X = num_training_set[:, :-1]
Y = num_training_set[:, -1]

X = torch.tensor(X, dtype=torch.float32)
Y = torch.tensor(Y, dtype=torch.float32).reshape(-1, 1)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=111)

mean = torch.mean(X_train, dim=0)
std = torch.std(X_train, dim=0)

X_train_normalized = (X_train - mean) / std
X_test_normalized = (X_test - mean) / std



model = nn.Sequential(
    nn.Linear(49, 74),
    nn.ReLU(),
    nn.Linear(74, 49),
    nn.ReLU(),
    nn.Linear(49, 1),
    nn.Sigmoid()
)



loss_fn = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)



# Use learning rate scheduler
scheduler = ExponentialLR(optimizer, gamma=0.95)

# training loop with learning rate scheduler and early stopping
best_loss = float('inf')
patience = 5
early_stop_counter = 0


# training

n_epochs = 100
batch_size = 32
 
for epoch in range(n_epochs):
    model.train()
    for i in range(0, len(X_train_normalized), batch_size):
        Xbatch = torch.tensor(X_train_normalized[i:i+batch_size], dtype=torch.float32)
        Ybatch = torch.tensor(Y_train[i:i+batch_size], dtype=torch.float32).reshape(-1, 1)
        Y_pred = model(Xbatch)
        loss = loss_fn(Y_pred, Ybatch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    # Update learning rate
    scheduler.step()
    
    # Validation
    model.eval()
    with torch.no_grad():
        Y_pred_val = model(torch.tensor(X_test_normalized, dtype=torch.float32))
        val_loss = loss_fn(Y_pred_val, torch.tensor(Y_test, dtype=torch.float32).reshape(-1, 1))
        
    print(f'Finished epoch {epoch}, training loss {loss}, validation loss {val_loss}')
    
    # Early stopping
    if val_loss < best_loss:
        best_loss = val_loss
        early_stop_counter = 0
    else:
        early_stop_counter += 1
        if early_stop_counter >= patience:
            print("Early stopping triggered.")
            break

# Evaluate model on test set
model.eval()
with torch.no_grad():
    Y_pred_test = model(torch.tensor(X_test_normalized, dtype=torch.float32))
    test_loss = loss_fn(Y_pred_test, torch.tensor(Y_test, dtype=torch.float32).reshape(-1, 1))
print(f"Test loss: {test_loss}")

model.state_dict()

torch.save(model.state_dict(), "neural_network_2403.pytorch")









## ROC

# Assuming y_test are 

# model = nn.Sequential(
#     nn.Linear(49, 74),
#     nn.ReLU(),
#     nn.Linear(74, 49),
#     nn.ReLU(),
#     nn.Linear(49, 1),
#     nn.Sigmoid()
# )

# model.load_state_dict(torch.load("neural_network_1303_2.pytorch"))

prediction = model(X_test_normalized)

# Compute the ROC curve
fpr, tpr, thresholds = roc_curve(Y_test, prediction.detach().numpy())

# Compute the AUC
roc_auc = auc(fpr, tpr)

# Plot the ROC curve
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()
# ...

# Evaluate model on test set
model.eval()
with torch.no_grad():
    Y_pred_test = model(torch.tensor(X_test_normalized, dtype=torch.float32))
    test_loss = loss_fn(Y_pred_test, torch.tensor(Y_test, dtype=torch.float32).reshape(-1, 1))
    Y_pred_test_binary = (Y_pred_test >= 0.5).float() # threshold ¿?¿?¿?¿
    report = classification_report(Y_test, Y_pred_test_binary)
    print(report)