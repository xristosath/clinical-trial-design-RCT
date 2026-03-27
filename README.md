# Clinical Trial Design in R

## 📌 Overview

This project implements key statistical methodologies used in **Phase I–II clinical trials**, including dose-escalation, binomial designs, randomization, and interim analysis.

The analysis is fully implemented in **R**, combining theoretical understanding with simulation-based validation.

---

## 🧪 Methods Implemented

### Phase I

* 3+3 Dose Escalation Design
* Transition Probabilities
* Maximum Tolerated Dose (MTD) estimation
* Target Dose Accuracy evaluation

### Phase II

* Exact Binomial Single-Stage Design
* Simon Two-Stage Design (Optimal & Minimax)
* Power analysis under intermediate efficacy

### Randomization

* Simple Randomization
* Block Randomization
* Imbalance analysis via simulation

### Interim Analysis

* O'Brien-Fleming boundaries
* Pocock boundaries
* Decision rules for early stopping

---

## 📊 Key Insights

* The **3+3 design** is conservative and tends to underestimate the MTD
* **Simon two-stage designs** reduce patient exposure under ineffective treatments
* **Block randomization** ensures balance but introduces predictability
* **O'Brien-Fleming** protects against early false positives, while **Pocock** allows earlier stopping

---

## 📁 Project Structure

```
R/         → R scripts
report/    → Full assignment report (PDF)
```

---

## ▶️ How to Run

1. Open the R script:

```
R/RCT_final.R
```

2. Install required packages:

```r
install.packages(c("clinfun", "blockrand", "gsDesign"))
```

3. Run the script step by step.

---

## 📄 Report

Full methodology, results, and interpretation are available in:

```
report/RCT_Final_Athanasakopoulos.pdf
```

---

## 🧠 Author

Christos Athanasakopoulos
MSc Biostatistics / Data Science

---
