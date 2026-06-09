

def mae_coverage_plotting():
    # Vollständiges Python-Skript zum Einlesen der Excel-Daten
    # und Erstellen der beiden Plots mit xticks nur an vorhandenen Sample Ratios

    import pandas as pd
    import matplotlib.pyplot as plt

    # Pfad zur Excel-Datei
    file_path = 'C:/Users/Frede/Nextcloud/Uni/Promotion/Projekte/OCEL Sampling/Results MAE coverage.xlsx'

    # Einlesen der Rohdaten ohne Header
    raw = pd.read_excel(file_path, sheet_name='rp mod', header=None)

    # Mapping der Algorithmus-Bezeichnungen
    algo_map = {
        'RP OCEL sampling': 'RP-OCEL',
        'random sampling': 'random sampling'
    }

    # Daten extrahieren
    data = []
    current_algo = None
    for _, row in raw.iterrows():
        algo_cell = row[1]
        if isinstance(algo_cell, str) and algo_cell in algo_map:
            current_algo = algo_map[algo_cell]
            continue
        if current_algo is None:
            continue
        try:
            sample_ratio = float(algo_cell)
            mae = float(row[2])
            coverage = float(row[3])
            new_metric = float(row[4])  # percentage of changed columns
            data.append({
                'Algorithmus': current_algo,
                'SampleRatio': sample_ratio,
                'MAE': mae,
                'Coverage': coverage,
                'NewMetric': new_metric
            })
        except (ValueError, TypeError):
            continue

    df = pd.DataFrame(data)
    xticks = sorted(df['SampleRatio'].dropna().unique())
    xticks = [str(x) for x in xticks]
    positions = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

    colors = {
        'RP-OCEL': 'red',
        'random sampling': 'blue'
    }

    # Plot 1: MAE vs. Sample Ratio
    plt.figure(figsize=(8, 5))
    for algo in ['RP-OCEL', 'random sampling']:
        subset = df[df['Algorithmus'] == algo]
        plt.plot(
            subset['SampleRatio'], subset['MAE'],
            marker='o', label=algo, color=colors[algo]
        )
    plt.xlabel('sample ratio', fontsize=15)
    plt.ylabel('MAE', fontsize=15)
    # plt.title('MAE vs. Sample Ratio für RP OCEL und Random Object')
    plt.xticks(positions, xticks, fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    plt.xlim(0.55, 0)
    plt.ylim(0, 320)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot 2: Coverage vs. Sample Ratio
    plt.figure(figsize=(8, 5))
    for algo in ['RP-OCEL', 'random sampling']:
        subset = df[df['Algorithmus'] == algo]
        plt.plot(
            subset['SampleRatio'], subset['Coverage'],
            marker='o', label=algo, color=colors[algo]
        )
    plt.xlabel('sample ratio', fontsize=15)
    plt.ylabel('coverage', fontsize=15)
    # plt.title('Coverage vs. Sample Ratio für RP OCEL und Random Object')
    plt.xticks(positions, xticks, fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    plt.xlim(0.55, 0)
    plt.ylim(-0.1, 1.1)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot 3: Percentage of Changed Columns vs. Sample Ratio
    plt.figure(figsize=(8.5, 5))
    for algo in ['RP-OCEL', 'random sampling']:
        subset = df[df['Algorithmus'] == algo]
        plt.plot(
            subset['SampleRatio'], subset['NewMetric'],
            marker='o', label=algo, color=colors[algo]
        )
    plt.xlabel('sample ratio', fontsize=15)
    plt.ylabel('percentage of events with \n original object relations (%)', fontsize=15)
    # plt.title('Percentage of Changed Columns vs. Sample Ratio für RP OCEL und Random Object')
    plt.xticks(positions, xticks, fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    plt.xlim( 0.55, 0)
    plt.ylim(-10, 110)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_fitness():
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load and clean data
    df = pd.read_excel('C:/Users/Frede/Downloads/FIrness_AB_vs_rd_self.xlsx', header=None)
    df = df.iloc[:, :3]
    df.columns = ['x', 'RP-OCEL', 'random sampling']
    df = df.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)

    positions = list(range(len(df)))
    # Replace 1.0 label with "reference"
    tick_labels = ["reference" if x == 1 else str(x) for x in df['x']]

    # # Plot 1: Raw values with y-axis 0-1, horizontal x-labels, larger fonts
    # plt.figure(figsize=(8, 5))
    # plt.plot(positions, df['RP-OCEL'], label='RP-OCEL', color='red', marker='o')
    # plt.plot(positions, df['random sampling'], label='Random sampling', color='blue', marker='o')
    # plt.ylim(-0.12, 1)
    # plt.xticks(positions, tick_labels, rotation=0, fontsize=18)
    # plt.yticks(fontsize=18)
    # plt.xlabel('sample ratio', fontsize=18)
    # plt.ylabel('fitness', fontsize=18)
    # # plt.title('AB-OCEL (red) and Random Sampling (blue) vs. x', fontsize=18)
    # plt.legend(fontsize=18)
    # plt.tight_layout()
    # plt.grid(True)
    # plt.show()

    # Compute reference values and deviations
    ref_ab = df.loc[df['x'] == 1, 'RP-OCEL'].iloc[0]
    ref_rs = df.loc[df['x'] == 1, 'random sampling'].iloc[0]
    abs_dev_ab = (df['RP-OCEL'] - ref_ab).abs()
    abs_dev_rs = (df['random sampling'] - ref_rs).abs()

    # Plot 2: Absolute deviations with horizontal x-labels, larger fonts
    plt.figure(figsize=(8, 5))
    plt.plot(positions, abs_dev_ab, label='RP-OCEL', color='red', marker='o')
    plt.plot(positions, abs_dev_rs, label='Random sampling', color='blue', marker='o')
    plt.ylim(-0.01, 0.5)
    plt.xticks(positions, tick_labels, rotation=0, fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('sample ratio', fontsize=15)
    plt.ylabel('deviation from reference', fontsize=15)
    # plt.title('Absolute Deviation (Reference at x=1)', fontsize=18)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.tight_layout()
    plt.show()





mae_coverage_plotting()
# plot_fitness()