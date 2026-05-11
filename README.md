# MSA-Survival-Prediction-Model
# MSA-Survival-Prediction-Model

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://dongxiao0502-msa-survial-streamapp-heudbv.streamlit.app/)

&gt; **Interpretable Machine Learning Model for Individualized Survival Prediction of Multiple System Atrophy**

This repository contains the implementation of a Random Survival Forest (RSF) model with SHAP-based interpretability for predicting survival outcomes in Multiple System Atrophy (MSA) patients.

---

## 🌐 Online Prediction Tool (Beta)

**Access the live application**: 
👉 **https://dongxiao0502-msa-survial-streamapp-heudbv.streamlit.app/**

### Quick Start Guide

| Step | Action | Details |
|:---|:---|:---|
| **1** | **Input Features** | Enter patient clinical features in the left sidebar |
| **2** | **Submit** | Click the "Submit" button |
| **3** | **View Results** | Review individualized predictions and explanations |

### Input Features (Left Sidebar)
- **Demographics**: Age
- **Clinical Symptoms**: Frequent Falls, Pathological Signs, Bradykinesia, Rigidity, Rest Tremor, Neurogenic OH, Urinary Incontinence, Constipation
- **Laboratory Values**: Residual Urine, Neutrophils, Lymphocytes, Uric Acid, Homocysteine, Vitamin B9

### Output Results
- 📈 **Predicted Survival Curve**: Visual survival probability over time
- ⏱️ **Median Survival Time**: Predicted survival duration (years)
- 📊 **Time-point Survival Rates**: 3-year, 5-year, 7-year survival probabilities
- 🔍 **SHAP-based Explanations**: Individualized feature contribution analysis

---

## 📁 Repository Structure

This repository is organized into two branches for different purposes:

### 🔬 Branch: `main` — Model Development & Optimization
**Purpose**: Model training, hyperparameter tuning, and comparative analysis

| File | Description | Language |
|:---|:---|:---|
| `survmodel0706.R` | Main RSF model implementation and training | R |
| `GBSA(1).txt` | GBSA (Gradient Boosting Survival Analysis) model parameters | Python |
| `RSF(1).txt` | Random Survival Forest model parameters | Python |
| `SVM(1).txt` | Support Vector Machine survival model parameters | Python |
| `coxnet(1).txt` | CoxNet (Cox proportional hazards with elastic net) parameters | Python |
| `xgb(1).txt` | XGBoost survival model parameters | Python |

**Note**: Original training data not uploaded due to privacy restrictions. Data available from corresponding author upon reasonable request.

---

### 🌐 Branch: `website` — Online Deployment
**Purpose**: Streamlit-based web application for clinical use

| File | Description |
|:---|:---|
| `streamapp.py` | Main Streamlit application code |
| `rsf_model.pkl` | Serialized trained RSF model |
| `feature_names.pkl` | Feature name mappings |
| `train111.csv` | Sample data structure template (synthetic/de-identified) |
| `requirements.txt` | Python dependencies for deployment |

**Switch to this branch**:
```bash
git checkout website
