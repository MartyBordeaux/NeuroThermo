from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(csv_path="component_boundary_summary_v1_3.csv", outdir="figures"):
    df = pd.read_csv(csv_path)
    b = df[(df["subset"] == "core_secure_pairs") & (df["projection"] == "isi")]
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    for stage in ["WT_exit", "balance", "SCA3_entry"]:
        fig, ax = plt.subplots(figsize=(8.5, 5.6))
        for mode in ["combined", "kappa_only", "J_only"]:
            # Explicit column indexing is required: DataFrame.mode is a pandas method.
            g = b[(b["stage"] == stage) & (b["mode"] == mode)].sort_values("p_intrinsic")
            ax.plot(
                g["p_intrinsic"],
                g["median_majority_support"],
                marker="o",
                markersize=3,
                linewidth=1.5,
                label=mode.replace("_", " "),
            )
        ax.set(
            xlabel="Intrinsic progress",
            ylabel="Required component progress",
            xlim=(0, 1),
            ylim=(0, 1),
            title=f"ISI {stage.replace('_', ' ')}: combined vs kappa-only vs J-only",
        )
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / f"01_boundary_{stage}.png", dpi=220)
        plt.close(fig)


if __name__ == "__main__":
    main()
