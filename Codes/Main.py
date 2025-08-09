
# ========================
# Main Script
# ========================

import torch
from torch_geometric.loader import DataLoader

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load your processed dataset
    train_dataset = torch.load("train_dataset.pt")
    val_dataset = torch.load("val_dataset.pt")
    test_dataset = torch.load("test_dataset.pt")

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    # Define model
    model = GNNModel(num_node_features=train_dataset.num_node_features, hidden_channels=64).to(device)

    # Train and evaluate
    model, y_true, y_pred = train_and_evaluate(model, train_loader, val_loader, test_loader, device)

    # Residual analysis
    residual_analysis(y_true, y_pred, save_path="residuals.png")

if __name__ == "__main__":
    main()
