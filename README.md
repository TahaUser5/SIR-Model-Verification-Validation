# SIR Epidemic Model — Verification & Validation

A Python-based SIR (Susceptible-Infected-Recovered) epidemic model simulator with formal property verification, exhaustive parameter tuning, and predictive validation against real-world Influenza data.

> **Course:** Introduction to Computational Modelling  
> **Institution:** University of the Punjab, College of Information Technology  
> **Instructor:** Dr. Syed Waqar ul Qounain Jaffry

---

## Overview

This program simulates the spread of infectious disease using the classic **SIR compartmental model**, solved numerically via the **Euler method**. It was applied to real Influenza spread data from the **USA** and **Netherlands** during 2009–2011 to demonstrate model verification, validation, and prediction.

### What It Does

| Feature | Description |
|---|---|
| **Trace Generation** | Simulate SIR dynamics with any custom parameters (N, I₀, R₀, β, γ) |
| **Property Verification** | Check 3 formal mathematical properties across 10,201 parameter combinations |
| **Parameter Tuning** | Exhaustive grid search over β and γ to find the best fit to training data |
| **Model Prediction** | Apply tuned parameters to unseen test data and report testing error |
| **Import / Export** | Export simulation traces to CSV; import empirical I(t) data from CSV |

---

## Mathematical Model

The SIR model divides a population of size **N** into three compartments:

```
dS/dt = -β * I * S / N
dI/dt =  β * I * S / N  -  γ * I
dR/dt =  γ * I
```

| Symbol | Meaning |
|---|---|
| S | Susceptible population |
| I | Infected population |
| R | Recovered population |
| β | Transmission rate (infection spread speed) |
| γ | Recovery rate (how fast people recover) |
| Δt | Euler integration time step (0.1) |

---

## Formal Properties Verified

Three logical properties are checked on every simulation trace:

| # | Property | Mathematical Definition |
|---|---|---|
| P1 | **Conservation** | S(t) + I(t) + R(t) = N for all t |
| P2 | **Non-negativity** | S(t) ≥ 0, I(t) ≥ 0, R(t) ≥ 0 for all t |
| P3 | **Monotonicity** | S(t+1) ≤ S(t) for all t |

---

## Dataset

Real-world weekly Influenza case counts from the USA and Netherlands:

| Dataset | Period | USA Population | NL Population | Weeks |
|---|---|---|---|---|
| Training | 2009–2010 | 60,000 | 5,000 | 55 |
| Testing | 2010–2011 | 60,000 | 5,000 | 55 |

### Parameter Search Space

| Parameter | Range | Step Size |
|---|---|---|
| β (transmission) | [0.00, 1.00] | 0.01 |
| γ (recovery) | [0.00, 1.00] | 0.01 |
| Δt (time step) | 0.1 | — |
| **Total combinations** | **10,201** | — |

---

## How to Run

### Prerequisites
- Python 3.6 or higher (no external libraries needed)

### Execution
```bash
python sir_model.py
```

### Menu Options
```
  1. Generate SIR Trace
  2. Check Formal Properties
  3. Parameter Tuning (Exhaustive Search)
  4. Model Prediction (Test Phase)
  5. Export / Import
  6. Exit
```

### Recommended Workflow
1. Run **Option 2** to verify formal properties
2. Run **Option 3** to find the best parameters via exhaustive search
3. Run **Option 4** to test the tuned model against 2010–2011 data

---

## Sample Results

### Best Parameters Found (Training Data 2009–2010)

| Country | β | γ | Training MSE |
|---|---|---|---|
| USA | 0.76 | 0.54 | 722,367.13 |
| Netherlands | 1.00 | 0.74 | 11,616.37 |

### Prediction Results (Test Data 2010–2011)

| Country | Testing MSE |
|---|---|
| USA | 1,968,936.01 |
| Netherlands | 3,166.06 |

---

## Generated Reports

Running the program automatically creates these report files:

| File | Contents |
|---|---|
| `Properties_Report.txt` | Formal property verification results |
| `Tuning_Report.txt` | Top 1% parameter settings with training errors |
| `Prediction_Report.txt` | Testing errors for the best-tuned model |

---

## Project Structure

```
├── sir_model.py           # Main program (single file, menu-driven)
├── Properties_Report.txt  # Auto-generated after Option 2
├── Tuning_Report.txt      # Auto-generated after Option 3
├── Prediction_Report.txt  # Auto-generated after Option 4
└── README.md              # This file
```

---

## Technologies

- **Language:** Python 3
- **Libraries:** `csv`, `sys` (standard library only — zero dependencies)
- **Numerical Method:** Euler forward integration

---

## License

This project is for educational purposes as part of CS-501 coursework.
