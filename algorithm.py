import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import torch.optim as optim
from sklearn.model_selection import train_test_split
from torch.optim.lr_scheduler import ExponentialLR
from sklearn.metrics import roc_curve, roc_auc_score, auc
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
from torchviz import make_dot


# Load the training set
training_set = pd.read_parquet('../trainingSet_20240326.parquet')


# Set "is_lbs" as the last column
training_set["is_lbs"] = training_set.pop("is_lbs")




# Convert the training set to a numpy array, and then to torch tensor

num_training_set = np.array(training_set)
X = num_training_set[:, :-1]
Y = num_training_set[:, -1]

X = torch.tensor(X, dtype=torch.float32)
Y = torch.tensor(Y, dtype=torch.float32).reshape(-1, 1)

# Split the data into training and test sets
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=111)

# Normalize the data using the mean and standard deviation of the training set

mean = torch.mean(X_train, dim=0)
std = torch.std(X_train, dim=0)

X_train_normalized = (X_train - mean) / std
X_test_normalized = (X_test - mean) / std




######################
### Neural network ###
######################



# Define the model
model = nn.Sequential(
    nn.Linear(50, 74),
    nn.ReLU(),
    nn.Linear(74, 50),
    nn.ReLU(),
    nn.Linear(50, 25),
    nn.ReLU(),
    nn.Linear(25, 1),
    nn.Sigmoid()
)

if __name__ == "__main__":
    # Define the loss function and optimizer

    loss_fn = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.005)



    # Use learning rate scheduler
    scheduler = ExponentialLR(optimizer, gamma=0.95)

    # Define early stopping parameters

    best_loss = float('inf')
    patience = 5
    early_stop_counter = 0


    # Training loop with early stopping

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

    # Save the model
    torch.save(model.state_dict(), "neural_network_2603_1988_pdbs_6.2A.pytorch")









    ## Load neural network model from state dict and assess model

    

    # model = nn.Sequential(
    #     nn.Linear(50, 74),
    #     nn.ReLU(),
    #     nn.Linear(74, 50),
    #     nn.ReLU(),
    #     nn.Linear(50, 25),
    #     nn.ReLU(),
    #     nn.Linear(25, 1),
    #     nn.Sigmoid()
    # )



    # model.load_state_dict(torch.load("neural_network_2503_1988_pdbs_martin.pytorch"))
    # model.state_dict()
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



    # Evaluate model on test set
    model.eval()
    with torch.no_grad():
        Y_pred_test = model(torch.tensor(X_test_normalized, dtype=torch.float32))
        test_loss = loss_fn(Y_pred_test, torch.tensor(Y_test, dtype=torch.float32).reshape(-1, 1))
        Y_pred_test_binary = (Y_pred_test >= 0.5).float() # threshold ¿?¿?¿?¿
        report = classification_report(Y_test, Y_pred_test_binary)
        print(report)


    # Get image of the model
        

    # Get image of the model
    output_to_plot = model(X_test_normalized)

    dot = make_dot(output_to_plot, params=dict(model.named_parameters()))
    dot.format = 'png'
    dot.render('../model')