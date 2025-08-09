
# ========================
# Train and Evaluation
# ========================

import torch
import torch.nn as nn
import torch.optim as optim
from torch_geometric.loader import DataLoader
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np

# Training function
def train(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    for data in loader:
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = criterion(output.view(-1), data.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)

# Evaluation function
def evaluate(model, loader, criterion, device):
    model.eval()
    total_loss = 0
    y_true, y_pred = [], []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            output = model(data)
            loss = criterion(output.view(-1), data.y)
            total_loss += loss.item()
            y_true.extend(data.y.cpu().numpy())
            y_pred.extend(output.view(-1).cpu().numpy())
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    r2 = r2_score(y_true, y_pred)
    return total_loss / len(loader), rmse, r2, np.array(y_true), np.array(y_pred)

# Main training loop
def train_and_evaluate(model, train_loader, val_loader, test_loader, device, lr=0.001, epochs=100):
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()
    best_val_rmse = float('inf')
    best_model_state = None

    for epoch in range(epochs):
        train_loss = train(model, train_loader, optimizer, criterion, device)
        val_loss, val_rmse, val_r2, _, _ = evaluate(model, val_loader, criterion, device)

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            best_model_state = model.state_dict()

        print(f"Epoch {epoch+1:03d} | "
              f"Train Loss: {train_loss:.4f} | "
              f"Val Loss: {val_loss:.4f} | "
              f"Val RMSE: {val_rmse:.4f} | "
              f"Val R²: {val_r2:.4f}")

    # Load best model
    model.load_state_dict(best_model_state)

    # Final test evaluation
    test_loss, test_rmse, test_r2, y_true, y_pred = evaluate(model, test_loader, criterion, device)
    print(f"Test Loss: {test_loss:.4f} | Test RMSE: {test_rmse:.4f} | Test R²: {test_r2:.4f}")

    return model, y_true, y_pred
