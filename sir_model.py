import csv
import sys

# ── Dataset: Influenza I(t) — Training 2009-2010 (55 weeks) ──
TRAIN_USA = [8,15,18,23,16,37,24,34,19,23,30,24,30,30,20,44,44,48,50,74,85,87,140,124,179,218,313,451,579,1011,1336,2155,3818,4733,4743,5186,4371,3545,2544,1943,1378,880,608,347,185,125,64,37,29,19,17,11,7,6,4]
TRAIN_NL  = [1,4,21,45,42,21,20,16,98,197,480,456,37,51,36,27,27,20,15,41,50,69,147,226,226,226,104,67,37,29,11,4,4,4,1,4,4,1,4,1,1,1,1,1,2,1,2,2,2,1,1,2,2,1,1]

# ── Dataset: Influenza I(t) — Testing 2010-2011 (55 weeks) ──
TEST_USA = [4,10,12,24,19,32,21,45,34,46,64,67,62,76,129,142,167,245,350,485,581,607,839,1433,1985,2355,2647,3022,3708,4424,4860,5024,5722,4553,3571,2657,1938,1340,909,627,391,192,125,62,35,18,16,11,5,9,2,4,4,4,7]
TEST_NL   = [1,1,1,1,1,2,1,1,2,1,1,1,1,1,1,2,3,1,2,3,5,5,8,24,32,53,153,174,160,86,122,79,52,32,21,13,11,6,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

# ── Model constants from assignment specification ──
USA_N = 60000
NL_N = 5000

BETA_RANGE  = [round(i * 0.01, 2) for i in range(101)]   # 0.00 to 1.00, step 0.01
GAMMA_RANGE = [round(i * 0.01, 2) for i in range(101)]   # 0.00 to 1.00, step 0.01
DT = 0.1                                                 # Euler integration step


# ── SIR Model (Euler Method) ──
def simulate(N, I0, R0, beta, gamma, weeks):
    S = N - I0 - R0
    Inf = float(I0)
    R = float(R0)
    
    trace = [{'week': 1, 'S': S, 'I': Inf, 'R': R}]
    steps = int(round(1.0 / DT))
    
    for w in range(2, weeks + 1):
        for _ in range(steps):
            dS = -beta * Inf * S / N * DT
            dI = (beta * Inf * S / N - gamma * Inf) * DT
            dR = gamma * Inf * DT
            
            S += dS
            Inf += dI
            R += dR
            
        trace.append({'week': w, 'S': S, 'I': Inf, 'R': R})
        
    return trace


# ── Three Formal / Logical Properties ──
def check_properties(trace, N):
    for i, pt in enumerate(trace):
        S = pt['S']
        Inf = pt['I']
        R = pt['R']
        
        # P1: Conservation - Total population remains constant
        if abs(S + Inf + R - N) > 1e-2:
            return False
            
        # P2: Non-negativity - No negative population counts
        if S < -1e-2 or Inf < -1e-2 or R < -1e-2:
            return False
            
        # P3: Monotonicity - Susceptible never increases
        if i > 0 and S > trace[i-1]['S'] + 1e-2:
            return False
            
    return True


# ── Helpers ──
def mse(trace, data):
    total_error = 0
    for i in range(len(data)):
        total_error += (trace[i]['I'] - data[i]) ** 2
    return total_error / len(data)


def export_csv(trace, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['week', 'S', 'I', 'R'])
        writer.writeheader()
        writer.writerows(trace)
    print(f"  [OK] Exported to '{filename}'")


def import_csv(filename):
    data = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append(float(row['I']))
    return data


# ── Exhaustive Search (shared by tuning options) ──
def exhaustive_search(data, N, I0, label):
    results = []
    print(f"\n  Searching for {label}...")
    
    for i, b in enumerate(BETA_RANGE):
        for g in GAMMA_RANGE:
            tr = simulate(N, I0, 0, b, g, len(data))
            error = mse(tr, data)
            results.append({'beta': b, 'gamma': g, 'error': error})
            
        sys.stdout.write(f"\r  Progress: {(i+1)/len(BETA_RANGE)*100:.1f}%")
        sys.stdout.flush()
        
    # Sort from smallest error to largest
    results.sort(key=lambda x: x['error'])
    return results


# ── Option 1: Generate a custom SIR trace ──
def opt_trace():
    print("\n" + "-"*50 + "\n  Option 1: Generate SIR Trace\n" + "-"*50)
    try:
        N = float(input("  Population (N): "))
        I0 = float(input("  Initial Infected I(0): "))
        R0 = float(input("  Initial Recovered R(0): "))
        beta = float(input("  Beta: "))
        gamma = float(input("  Gamma: "))
        weeks = int(input("  Weeks: "))
    except ValueError:
        print("  [ERROR] Enter numeric values.")
        return

    trace = simulate(N, I0, R0, beta, gamma, weeks)

    print(f"\n  {'Week':<7} {'S':<16} {'I':<16} {'R':<16}")
    print("  " + "-"*55)
    for t in trace:
        print(f"  {t['week']:<7} {t['S']:<16.4f} {t['I']:<16.4f} {t['R']:<16.4f}")

    if input("\n  Export to CSV? (y/n): ").lower() == 'y':
        filename = input("  Filename: ")
        export_csv(trace, filename)


# ── Option 2: Check formal properties on all parameter sets ──
def opt_properties():
    print("\n" + "-"*50 + "\n  Option 2: Check Formal Properties\n" + "-"*50)
    print("  P1: Conservation  — S+I+R = N")
    print("  P2: Non-negativity — S,I,R >= 0")
    print("  P3: Monotonicity  — S never increases")

    total = len(BETA_RANGE) * len(GAMMA_RANGE)
    weeks = len(TRAIN_USA)
    counts = {}

    test_configurations = [
        ("USA", USA_N, TRAIN_USA[0]), 
        ("NL", NL_N, TRAIN_NL[0])
    ]

    for label, N, I0 in test_configurations:
        passed = 0
        print(f"\n  Checking {label}...")
        
        for i, b in enumerate(BETA_RANGE):
            for g in GAMMA_RANGE:
                trace = simulate(N, I0, 0, b, g, weeks)
                if check_properties(trace, N):
                    passed += 1
                    
            sys.stdout.write(f"\r  Progress: {(i+1)/len(BETA_RANGE)*100:.1f}%")
            sys.stdout.flush()
            
        counts[label] = (passed, total - passed)

    print(f"\n\n  Results ({total} parameter sets each):")
    for label, (p, f_) in counts.items():
        print(f"  {label:15}: {p} passed, {f_} failed")

    with open("Properties_Report.txt", "w") as f:
        f.write("FORMAL PROPERTIES VERIFICATION REPORT\n" + "="*50 + "\n\n")
        f.write("P1 (Conservation) : S + I + R = N at all times\n")
        f.write("P2 (Non-negativity): S, I, R >= 0 at all times\n")
        f.write("P3 (Monotonicity) : S(t+1) <= S(t) for all t\n\n")
        f.write(f"Parameter Space: Beta [0.00,1.00] step 0.01, Gamma [0.00,1.00] step 0.01\n")
        f.write(f"Total Combinations: {total}\n\n")
        for label, (p, fl) in counts.items():
            f.write(f"{label}: Passed {p}/{total}, Failed {fl}/{total}\n")
            
    print("  [OK] Report saved to 'Properties_Report.txt'")


# ── Option 3: Parameter tuning via exhaustive search ──
def opt_tuning():
    print("\n" + "-"*50 + "\n  Option 3: Parameter Tuning\n" + "-"*50)
    print("  1. Use built-in training data (2009-2010)")
    print("  2. Use imported CSV data")

    choice = input("  Choice (1/2): ")
    if choice == '2':
        N = float(input("  Population (N): "))
        I0 = float(input("  Initial Infected I(0): "))
        filename = input("  CSV filename (columns: week,I): ")
        data = import_csv(filename)
        
        res = exhaustive_search(data, N, I0, "Imported")
        print(f"\n  Best: Beta={res[0]['beta']}, Gamma={res[0]['gamma']} (MSE={res[0]['error']:.2f})")
        return res[0], None

    usa_res = exhaustive_search(TRAIN_USA, USA_N, TRAIN_USA[0], "USA")
    nl_res = exhaustive_search(TRAIN_NL, NL_N, TRAIN_NL[0], "Netherlands")
    top = max(1, int(len(usa_res) * 0.01))

    test_results = [
        ("USA", usa_res), 
        ("Netherlands", nl_res)
    ]

    with open("Tuning_Report.txt", "w") as f:
        f.write("PARAMETER TUNING REPORT — Top 1%\n" + "="*50 + "\n\n")
        for label, res in test_results:
            f.write(f"{label} — Top 1% ({top} of {len(res)}):\n")
            f.write(f"{'Rank':<6}| {'Beta':<8}| {'Gamma':<8}| {'MSE':<15}\n" + "-"*40 + "\n")
            for i in range(top):
                r = res[i]
                f.write(f"{i+1:<6}| {r['beta']:<8.2f}| {r['gamma']:<8.2f}| {r['error']:<15.4f}\n")
            f.write("\n")

    print(f"\n  USA Best : Beta={usa_res[0]['beta']}, Gamma={usa_res[0]['gamma']} (MSE={usa_res[0]['error']:.2f})")
    print(f"  NL Best  : Beta={nl_res[0]['beta']}, Gamma={nl_res[0]['gamma']} (MSE={nl_res[0]['error']:.2f})")
    print("  [OK] Report saved to 'Tuning_Report.txt'")
    
    return usa_res[0], nl_res[0]


# ── Option 4: Model prediction on test data ──
def opt_prediction(best_usa, best_nl):
    print("\n" + "-"*50 + "\n  Option 4: Model Prediction\n" + "-"*50)
    if not best_usa or not best_nl:
        print("  [!] Run Option 3 first.")
        return

    results = {}
    test_configurations = [
        ("USA", USA_N, TEST_USA[0], TEST_USA, best_usa),
        ("NL", NL_N, TEST_NL[0], TEST_NL, best_nl)
    ]
    
    for label, N, I0, test, best in test_configurations:
        tr = simulate(N, I0, 0, best['beta'], best['gamma'], len(test))
        results[label] = {
            'trace': tr, 
            'mse': mse(tr, test), 
            'best': best
        }

    print("\n  Prediction Results (Test Data 2010-2011):")
    for label, r in results.items():
        print(f"  {label:5}: MSE = {r['mse']:.4f}  (Beta={r['best']['beta']}, Gamma={r['best']['gamma']})")

    with open("Prediction_Report.txt", "w") as f:
        f.write("MODEL PREDICTION REPORT — Test Data 2010-2011\n" + "="*50 + "\n\n")
        for label, r in results.items():
            f.write(f"{label}:\n  Beta={r['best']['beta']}, Gamma={r['best']['gamma']}\n  Testing MSE = {r['mse']:.4f}\n\n")
            
    print("  [OK] Report saved to 'Prediction_Report.txt'")

    if input("\n  Export prediction traces to CSV? (y/n): ").lower() == 'y':
        for label, r in results.items():
            filename = f"prediction_{label.lower()}.csv"
            export_csv(r['trace'], filename)


# ── Option 5: Import / Export utilities ──
def opt_io():
    print("\n" + "-"*50 + "\n  Option 5: Export / Import\n" + "-"*50)
    print("  1. Export a SIR trace to CSV")
    print("  2. Import empirical I(t) from CSV")
    print("  3. Back")
    
    c = input("  Choice (1/2/3): ")
    
    if c == '1':
        try:
            N = float(input("  Population (N): "))
            I0 = float(input("  I(0): "))
            R0 = float(input("  R(0): "))
            beta = float(input("  Beta: "))
            gamma = float(input("  Gamma: "))
            weeks = int(input("  Weeks: "))
            fn = input("  Output filename: ")
            
            trace = simulate(N, I0, R0, beta, gamma, weeks)
            export_csv(trace, fn)
        except ValueError:
            print("  [ERROR] Invalid input.")
            
    elif c == '2':
        filename = input("  CSV filename (columns: week,I): ")
        data = import_csv(filename)
        print(f"  Imported {len(data)} weeks. Preview: {data[:5]}")


# ── Main Menu ──
def main():
    best_usa = None
    best_nl = None
    
    while True:
        print("\n  1. Generate SIR Trace")
        print("  2. Check Formal Properties")
        print("  3. Parameter Tuning (Exhaustive Search)")
        print("  4. Model Prediction (Test Phase)")
        print("  5. Export / Import")
        print("  6. Exit")
        
        c = input("\n  Choice (1-6): ")
        
        if c == '1': 
            opt_trace()
        elif c == '2': 
            opt_properties()
        elif c == '3': 
            best_usa, best_nl = opt_tuning()
        elif c == '4': 
            opt_prediction(best_usa, best_nl)
        elif c == '5': 
            opt_io()
        elif c == '6': 
            print("  Exit")
            break
        else: 
            print("  [ERROR] Invalid choice.")

if __name__ == "__main__":
    main()
