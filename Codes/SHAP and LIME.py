
# ========================
# SHAP Explainability and Mean Bar Plot
# ========================

import shap
import numpy as np
import matplotlib.pyplot as plt

 # ---------------------------
 # SHAP Explainability Section - On All Test Data
 # ---------------------------
 print("Running SHAP analysis on all train data...")
 X_shap = extract_graph_features(train_loader)
 X_tensor = torch.tensor(X_shap, dtype=torch.float)

 mlp_model = model.post_pool.cpu()
 def model_predict(x_numpy):
       x_tensor = torch.tensor(x_numpy, dtype=torch.float)
       with torch.no_grad():
           preds = mlp_model(x_tensor)
       return preds.cpu().numpy().flatten()

 feature_names = [f"GNN_feat_{i}" for i in range(64)] + ["MolWt", "TPSA", "NumRings", "LogP", "ValenceElectrons"]

 explainer = shap.KernelExplainer(model_predict, X_shap[:100])  # Background subset for efficiency
 shap_values = explainer.shap_values(X_shap)  # Full test set

 shap.summary_plot(shap_values, X_shap, feature_names=feature_names, show=False, sort=True)
 plt.savefig(os.path.join(output_dir, "shap_summary_plot_full_ALL.png"), dpi=300)
 plt.close()

 print("SHAP analysis completed on all test data.")
 
 
   # ---------------------------
 # SHAP Mean Bar Plot
 # ---------------------------
 print("Plotting SHAP mean feature importance...")
 shap_mean = np.mean(shap_values, axis=0)
 plt.figure(figsize=(10, 6))
 plt.barh(feature_names, shap_mean, color='cornflowerblue')
 plt.xlabel("Mean SHAP Value")
 plt.title("Mean SHAP Value across Test Set")
 plt.tight_layout()
 plt.savefig(os.path.join(output_dir, "shap_mean_bar_plot.png"), dpi=300)
 plt.close()


   # ---------------------------
   # LIME on 50 samples in loop
   # ---------------------------
   from lime.lime_tabular import LimeTabularExplainer
   
     print("Running LIME analysis on 50 samples...")
     explainer_lime = lime_tabular.LimeTabularExplainer(
         training_data=X_shap,
         mode='regression',
         feature_names=feature_names,
         discretize_continuous=True
     )

     for i in range(50):
         explanation = explainer_lime.explain_instance(
             X_shap[i],
             model_predict,
             num_features=10
         )
         fig = explanation.as_pyplot_figure()
         fig.tight_layout()
         fig.savefig(os.path.join(output_dir, f"lime_explanation_sample{i}.png"), dpi=300)
         plt.close(fig)
         explanation.save_to_file(os.path.join(output_dir, f"lime_explanation_sample{i}.html"))

     print("LIME analysis completed on 50 samples.")
