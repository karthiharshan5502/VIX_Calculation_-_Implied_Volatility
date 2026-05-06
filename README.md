# VIX_Calculation_-_Implied_Volatility
Implementing the CBOE VIX formula and at-the-money implied volatility estimation using S&amp;P 500 options data in SAS (Jan–Feb 2015).


## 📌 Overview

This project replicates the **CBOE Volatility Index (VIX)** methodology from scratch using real S&P 500 options data and estimates **at-the-money implied volatility** via the **Newton-Raphson method**. Implemented entirely in **SAS (IML)** as part of MATH802 — Advanced Financial Modelling & Analytics at AUT (Semester 1, 2024).

---

## 🔬 Project Structure

The project is divided into three parts:

**Part 1 — VIX Calculation**
Implement the CBOE VIX formula using near-term and next-term S&P 500 option data for a chosen date, following the official CBOE white paper methodology.

**Part 2 — At-the-Money Implied Volatility**
Estimate implied volatility at the strike price closest to the market price using Newton-Raphson iteration applied to the Black-Scholes model.

**Part 3 — Comparative Analysis**
Calculate VIX and implied volatility across multiple dates (Jan–Feb 2015) and write an analysis comparing their differences and similarities.

---

## 📂 Repository Structure

```
├── MATH802_project_2024.pdf        # Full written project report with SAS outputs
├── MATH802.pdf                     # Assignment brief and instructions
└── vixwhite.pdf                    # CBOE VIX White Paper (reference methodology)
```

---

## 📊 Dataset

- **Source**: `option(Jan2015-Feb2015).csv` — S&P 500 options data (provided via course Canvas)
- **Date range**: January 2015 – February 2015
- **Columns**:

| Column | Description |
|---|---|
| `date` | Date of the options |
| `tau` | Time to option maturity (calendar days) |
| `K` | Option strike price |
| `call_option_price` | Call option price |
| `S` | S&P 500 index value on the day |
| `r` | Annual risk-free interest rate |
| `put_option_price` | Put option price |

---

## ⚙️ Methodology

### Part 1 — VIX Calculation (CBOE Formula)

**Step 1** — Compute time to expiration T₁ (near-term) and T₂ (next-term):
```
T = tau / 365
T1 = 15/365 = 0.0410959   (near-term)
T2 = 35/365 = 0.0958904   (next-term)
```

**Step 2** — Calculate forward index prices F₁ and F₂:
```
F = K + e^(R×T) × (Call Price − Put Price)
```
Chosen at the strike where |Call − Put| is minimised.

**Step 3** — Determine at-the-money strike K₀ (largest K ≤ F).

**Step 4** — Calculate mid-quote prices:
```
Midquote = (Call Price + Put Price) / 2
```

**Step 5** — Compute ΔKᵢ (interval between strikes) and contribution by strike:
```
Contribution = (ΔKᵢ / Kᵢ²) × e^(R×T) × Q(Kᵢ)
```

**Step 6** — Calculate near- and next-term variance (σ²):
```
σ² = (2/T) × Σ Contribution − (1/T) × [(F/K₀) − 1]²
```

**Step 7** — Compute 30-day weighted average and final VIX:
```
VIX = 100 × √(WA)
```

### Part 2 — Implied Volatility (Newton-Raphson)

Using Black-Scholes for call price and iteratively solving for σ:
```
σ_new = σ - (C_BS - C_market) / Vega
```
Convergence achieved in 3 iterations.

---

## 📈 Key Results

| Date | VIX | Near-Term IV | Next-Term IV |
|---|---|---|---|
| 02 Jan 2015 | 20.09 | 0.1302 | 0.1369 |
| 15 Jan 2015 | 24.99 | 0.2158 | 0.1775 |
| 30 Jan 2015 | 25.19 | 0.1883 | 0.1703 |
| 02 Feb 2015 | 22.51 | 0.1447 | 0.1488 |
| 17 Feb 2015 | 19.72 | 0.1079 | 0.1203 |
| 26 Feb 2015 | 18.53 | 0.0929 | 0.1100 |

**Key findings:**
- VIX and both near/next-term IV trend together — higher VIX signals higher expected volatility
- Near-term IV reacts more sharply to immediate market conditions than next-term IV
- At low VIX levels, next-term IV exceeds near-term IV; this reverses at higher VIX levels

---

## 🛠️ Requirements

- **SAS** (with IML module)
- S&P 500 options data CSV (`option(Jan2015-Feb2015).csv`)

---

## 🚀 How to Run

```sas
/* 1. Import and filter data for a specific date */
proc import out=option
  datafile="/path/to/option(Jan2015-Feb2015).csv"
  dbms="csv";
  delimiter=',';
  getnames='yes';
run;

/* 2. Run the PROC IML block to calculate VIX */
/* 3. Run Newton-Raphson IML block for implied volatility */
/* Full code available in MATH802_project_2024.pdf */
```

The SAS code is written generically — simply change the date filter to calculate VIX for any available date in the dataset.

---

## 📄 References

- CBOE (2014). *The CBOE Volatility Index — VIX White Paper*. Chicago Board Options Exchange.
- Black, F. & Scholes, M. (1973). The pricing of options and corporate liabilities. *Journal of Political Economy*, 81(3), 637–654.

---

## 🏫 Course Information

| Field | Detail |
|---|---|
| Course | MATH802 — Advanced Financial Modelling & Analytics |
| Institution | Auckland University of Technology (AUT) |
| Semester | Semester 1, 2024 |
| Assessment | Written Project (40%) |
