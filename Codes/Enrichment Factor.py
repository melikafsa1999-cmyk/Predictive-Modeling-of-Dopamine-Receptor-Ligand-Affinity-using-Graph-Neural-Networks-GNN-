
# ========================
# Enrichment Factor (EF)
# ========================

import numpy as np

 # ---------------------------
 # External SMILES Prediction on full shuffled file
 # ---------------------------
 print("Running external prediction on entire dataset...")
 
 full_dataset = SMILESDataset(df) 
 full_loader = DataLoader(full_dataset, batch_size=32, shuffle=False)
 model.eval()
 full_preds = []
 
 with torch.no_grad():
     for data in full_loader:
         data = data.to(device)
         out = model(data.x, data.edge_index, data.edge_attr, data.batch, data.u)
         full_preds.append(out.cpu().numpy())
 
 full_preds = np.concatenate(full_preds, axis=0)
 df['Predicted'] = full_preds
 df.to_excel(os.path.join(output_dir, "shuffled_df_with_predictions.xlsx"), index=False)
 
 print("Prediction completed. Saved to 'shuffled_df_with_predictions.xlsx'.")
 
     
 # ---------------------------
 # Enrichment Factor (EF) Calculation
 # ---------------------------
 print("Calculating Enrichment Factor (EF)...")
 
 df_sorted = df.copy()
 df_sorted = df_sorted.sort_values(by="Predicted", ascending=False).reset_index(drop=True)
 
 top_k_percent = 0.1  
 k = max(1, int(top_k_percent * len(df_sorted)))
 
 # Top-k predicted ligands
 top_k = df_sorted.iloc[:k]
 

 affinity_threshold = df_sorted["pAffinity"].quantile(0.9)
 df_sorted["IsActive"] = (df_sorted["pAffinity"] >= affinity_threshold).astype(int)
 top_k["IsActive"] = (top_k["pAffinity"] >= affinity_threshold).astype(int)
 
 num_actives_total = df_sorted["IsActive"].sum()
 num_inactives_total = len(df_sorted) - num_actives_total
 num_actives_topk = top_k["IsActive"].sum()
 num_inactives_topk = k - num_actives_topk
 
 EF_k = (num_actives_topk / k) / (num_actives_total / len(df_sorted))
 
 print(f"Enrichment Factor @ top {int(top_k_percent * 100)}%: {EF_k:.3f}")
 

 ef_output_path = os.path.join(output_dir, f"enrichment_factor_top_{int(top_k_percent * 100)}.txt")
 with open(ef_output_path, "w") as f:
     f.write(f"Enrichment Factor @ top {int(top_k_percent * 100)}%: {EF_k:.3f}\n")
     f.write(f"Top {k} samples (out of {len(df_sorted)} total)\n")
     f.write(f"Total actives: {num_actives_total}\n")
     f.write(f"Actives in top {k}: {num_actives_topk}\n")
     f.write(f"Total inactives: {num_inactives_total}\n")
     f.write(f"Inactives in top {k}: {num_inactives_topk}\n")
 
 print(f"EF written to {ef_output_path}")