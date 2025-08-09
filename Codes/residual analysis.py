
# ========================
# Residual Analysis
# ========================
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_residuals_vs_predicted(preds, trues, idx, output_dir):
    residuals = trues - preds
    std_residuals = np.std(residuals)
    plt.figure(figsize=(6, 5))
    plt.scatter(preds, residuals, alpha=0.7)
    for i in range(len(preds)):
        if residuals[i]> 3*std_residuals or residuals[i]<-3*std_residuals:
            plt.annotate(idx[i], (preds[i], residuals[i]), textcoords="offset points", xytext=(5,5), ha='center', fontsize = 12)
    plt.axhline(y=0, color='red', linestyle='--')
    plt.xlabel('Predicted')
    plt.ylabel('Residuals')
    plt.title('Residuals vs. Predicted')
    plt.savefig(os.path.join(output_dir, "Residuals_vs_Predicted.png"), dpi=300)
    plt.show()
    
def residuals_histogram(preds, trues, output_dir, dataset_label):
    """Plot and save histogram of residuals with frequency counts"""
    residuals = trues - preds
    
    plt.figure(figsize=(8, 6))
    sns.histplot(
        residuals,
        kde=True,
        bins=30,
        stat='count',  # Frequency counts
        color='blue',
        edgecolor='black',
        alpha=0.7
    )
    
    # Add distribution parameters to title
    mu, std = np.mean(residuals), np.std(residuals)
    plt.title(f'Histogram of Residuals - {dataset_label} Set\nμ={mu:.2f}, σ={std:.2f}')
    plt.xlabel('Residuals')
    plt.ylabel('Frequency')
    
    # Save with unique filename
    filename = os.path.join(output_dir, f"residuals_hist_{dataset_label}.png")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Prevents plot from displaying twice if in notebook

def normal_probability_plot(preds, trues, output_dir):
    residuals = trues - preds
    plt.figure(figsize=(6,5))
    stats.probplot(residuals, dist="norm", plot=plt)
    plt.title('Normal Probability Plot of Residuals')
    plt.savefig(os.path.join(output_dir, "NPP.png"), dpi=300)
    plt.show()

def compute_hat_matrix_approx(features,preds, trues, idx):
    """
    Approximate the Hat matrix diagonal by using features in a linear sense.
    features: N x D (where N is the number of samples, D is dimension).
    Returns the approximate leverage for each sample (diagonal of the Hat matrix).
    """
    # In conventional linear regression, H = X (X'X)^-1 X'.
    # We'll do that with 'features' here, though GNN usage is more complex.
    X = np.array(features)
    n_vars=X.shape[0]
    n_samples=X.shape[1]
    # Add a bias column if you want to incorporate an intercept in the design
    # Example:
    # X = np.hstack([X, np.ones((X.shape[0], 1))])
    XtX = X.T @ X
    try:
        XtX_inv = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        # In case of singular matrix, use pseudo-inverse
        XtX_inv = np.linalg.pinv(XtX)
    H = X @ XtX_inv @ X.T
    leverages = np.diag(H)
    # The diagonal of H is the leverage for each sample
    critical_leverages =  3 * (n_vars + 1) / n_samples
    
    residuals = trues - preds
    std_resid = (residuals - np.mean(residuals)) / (np.std(residuals) + 1e-9)
    plt.figure(figsize=(6,5))
    plt.scatter(leverages, std_resid, alpha=0.7)

    plt.figure()
    plt.scatter(leverages, std_resid, alpha=0.6)
    plt.axhline(y=3, color='r', linestyle='--', label='±3 Std Dev Limit')
    plt.axhline(y=-3, color='r', linestyle='--')
    plt.axvline(x=critical_leverages, color='g', linestyle='--', label='Leverage Limit')
    for i in range(len(leverages)):
        plt.annotate(idx[i], (leverages[i], std_resid[i]),textcoords= "offset points", xytext= (5,5), ha='center', fontsize=12)
    plt.xlabel('Leverage')
    plt.ylabel('Standardized Residuals')
    plt.title('Williams Plot')
    plt.legend()
    plt.grid(True)
    plt.show()

    return 
    

def williams_plot(train, val, test, train_pred, val_pred, test_pred, train_nom, val_nom, test_nom, train_idx, val_idx, test_idx, output_dir):
    """
    plot all three data sets on one williams plot
    """
    #Notation
    #def notation(h, h_crit, std_residuals, idx):
        #for i in range(len(h)):
            #if h[i]>h_crit or std_residuals[i]>3 or std_residuals[i]<-3:
                #plt.annotate(idx[i], (h[i],std_residuals[i]),textcoords= "offset points", xytext= (5,5), ha='center', fontsize=12)
              
    #Calculate leverage values
    try:
        scores_dict = {
            'train' : train,
            'val' : val,
            'test' : test
            }
        train_inv = np.linalg.pinv(train.T @ train)
        leverage_dict = {}
        for dset, scores in scores_dict.items():
            leverage_dict[dset] = np.diag(scores @ train_inv @ scores.T)
              
        train_lev = leverage_dict['train']
        val_lev = leverage_dict['val']
        test_lev = leverage_dict['test']
        
    except Exception as e:
        print(f'the error in leverage calculation is: {e}')
    
    #Calculate std residuals
    try:
        pred_dict = {
            'train':{'pred': train_pred, 'nom': train_nom},
            'val':{'pred':val_pred, 'nom': val_nom},
            'test':{'pred':test_pred, 'nom': test_nom}
            }

        std_residual_dict = {}  
        for name, data in pred_dict.items():
            residual = data['nom'] - data['pred']
            std_residual_dict[name] = (residual - np.mean(residual))/(np.std(residual) + 1e-9)
            
        train_std_residuals = std_residual_dict['train']
        val_std_residuals = std_residual_dict['val']
        test_std_residuals = std_residual_dict['test']
    except Exception as e:
        print(f"The error in std_residual calculation is: {e}")
        
    # Critical leverage threshold
    p = train.shape[1]  # Number of variables
    n = len(train_pred)
    h_crit = 3 * (p + 1) / n
    
    plt.figure(figsize=(10,8))
    plt.scatter(train_lev, train_std_residuals, label= 'Train', edgecolors='k', color = 'blue', s = 60)
    plt.scatter(val_lev, val_std_residuals, label= 'Validation', edgecolors= 'k', color = 'green', s = 60)
    plt.scatter(test_lev, test_std_residuals, label = 'Test', marker= 'D',edgecolors='k', color = 'orange', s = 60)
    #notation(train_lev, h_crit, train_std_residuals, train_idx)
    #notation(val_lev, h_crit, val_std_residuals, val_idx)
    #notation(test_lev, h_crit, test_std_residuals, test_idx)
    plt.axhline(y=3, color='r' , linestyle = '--')
    plt.axhline(y=-3, color= 'r', linestyle = '--')
    plt.axvline(x = h_crit, color= 'g', linestyle = '--')
    plt.annotate(f' h* = {h_crit:0.2f}', (h_crit, np.min(train_std_residuals)), fontsize = 16)
    plt.legend(fontsize = 14, loc = 'upper right')
    plt.xlabel('Leverage', fontsize = 14)
    plt.ylabel('std residuals', fontsize = 14)
    plt.title('Williams plot', fontsize = 16)
    plt.savefig(os.path.join(output_dir,'Williams_Plot.png'), dpi=300)
    plt.show()
    
def evaluate_statistics_for_datasets(model, train_loader, val_loader, test_loader, device):
    """
    Computes typical regression statistics (MSE, R², slope, intercept) 
    for three datasets: train, val, test.
    Returns a dictionary of results.
    """
    results = {}
    for name, loader in zip(['Train', 'Val', 'Test'], [train_loader, val_loader, test_loader]):
        avg_loss, r2, slope, intercept, preds, trues = evaluate(model, loader, device)
        mae = mean_absolute_error(trues, preds)
        results[name] = {
            'MSE': avg_loss,
            'MAE': mae,
            'R2': r2,
            'Slope': slope,
            'Intercept': intercept
        }
    return results