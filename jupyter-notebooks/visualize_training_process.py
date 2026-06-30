#!/usr/bin/env python3
"""
Visualize PyTorch Lightning training process from CSV log files.

Automatically detects metric granularity:
- Epoch-level metrics (one value per epoch, e.g., val_loss) → x-axis: epoch
- Step-level metrics (multiple values per epoch, e.g., train_loss_step) → x-axis: step

Usage:
    python visualize_training_process.py <csv_file> [csv_file2, ...]
    python visualize_training_process.py --help
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def load_csv_files(csv_paths: List[str]) -> List[Tuple[str, pd.DataFrame]]:
    """Load one or more CSV files and return with labels."""
    data = []
    for path in csv_paths:
        csv_path = Path(path)
        if not csv_path.exists():
            print(f"Warning: {path} not found, skipping")
            continue

        df = pd.read_csv(csv_path)
        label = csv_path.stem  # Use filename without extension as label
        data.append((label, df))
        print(f"Loaded {path}: {len(df)} rows, {len(df.columns)} columns")

    return data


def detect_metric_columns(df: pd.DataFrame) -> dict:
    """
    Detect and categorize metric columns from PyTorch Lightning CSV.
    Also detect granularity (epoch-level or step-level) for each metric.

    Returns:
        dict: Categorized columns with keys 'loss', 'accuracy', 'learning_rate', 'other'
              Each metric entry includes granularity info
    """
    columns = df.columns.tolist()
    metrics = {
        'loss': [],
        'accuracy': [],
        'learning_rate': [],
        'other': []
    }

    # Common PyTorch Lightning metric patterns
    loss_patterns = ['loss', 'Loss']
    acc_patterns = ['acc', 'accuracy', 'Accuracy']
    lr_patterns = ['lr', 'learning_rate', 'learning-rate']

    for col in columns:
        if col in ['epoch', 'step']:
            continue

        # Detect granularity
        granularity = detect_granularity(df, col)

        metric_info = {'name': col, 'granularity': granularity}

        # Categorize based on patterns
        if any(pattern in col.lower() for pattern in loss_patterns):
            metrics['loss'].append(metric_info)
        elif any(pattern in col.lower() for pattern in acc_patterns):
            metrics['accuracy'].append(metric_info)
        elif any(pattern in col.lower() for pattern in lr_patterns):
            metrics['learning_rate'].append(metric_info)
        else:
            metrics['other'].append(metric_info)

    return metrics


def detect_granularity(df: pd.DataFrame, metric_col: str) -> str:
    """
    Detect if a metric is epoch-level or step-level.

    Epoch-level: one value per epoch (e.g., val_loss, val_accuracy)
    Step-level: multiple values per epoch (e.g., train_loss_step)

    Returns:
        str: 'epoch' or 'step'
    """
    if 'epoch' not in df.columns or metric_col not in df.columns:
        return 'step'  # Default to step if no epoch column

    # Get non-NaN values for this metric
    metric_df = df[['epoch', metric_col]].dropna()

    if len(metric_df) == 0:
        return 'step'

    # Count values per epoch
    values_per_epoch = metric_df.groupby('epoch').size()

    # If most epochs have exactly 1 value, it's epoch-level
    # Otherwise it's step-level
    if values_per_epoch.mode().iloc[0] == 1 and (values_per_epoch == 1).sum() / len(values_per_epoch) > 0.8:
        return 'epoch'
    else:
        return 'step'


def plot_single_metric(data: List[Tuple[str, pd.DataFrame]],
                       metric_info: dict,
                       ax: plt.Axes,
                       show_legend: bool = True):
    """
    Plot a single metric across multiple runs.

    Args:
        data: List of (label, dataframe) tuples
        metric_info: dict with 'name' and 'granularity' keys
        ax: matplotlib axes to plot on
        show_legend: whether to show legend
    """
    metric_name = metric_info['name']
    granularity = metric_info['granularity']

    # Choose x-axis based on granularity
    x_axis = 'epoch' if granularity == 'epoch' else 'step'

    for label, df in data:
        if metric_name not in df.columns:
            continue

        if x_axis not in df.columns:
            print(f"Warning: {x_axis} column not found in {label}")
            continue

        # Drop NaN values for this metric
        plot_df = df[[x_axis, metric_name]].dropna()

        if len(plot_df) > 0:
            ax.plot(plot_df[x_axis], plot_df[metric_name],
                   marker='o', markersize=3, label=label, linewidth=1.5)

    ax.set_xlabel(x_axis.capitalize())
    ax.set_ylabel(metric_name)
    ax.set_title(f"{metric_name} ({granularity}-level)")
    ax.grid(True, alpha=0.3)

    if show_legend and len(data) > 1:
        ax.legend()


def plot_training_curves(data: List[Tuple[str, pd.DataFrame]],
                         output_path: Optional[str] = None,
                         figsize: Tuple[int, int] = (15, 10)):
    """
    Create comprehensive training visualization.

    Args:
        data: List of (label, dataframe) tuples
        output_path: Path to save the figure
        figsize: Figure size (width, height)
    """
    if not data:
        print("No data to plot")
        return

    # Detect metrics from the first dataframe
    metrics = detect_metric_columns(data[0][1])
    print(f"\nDetected metrics:")
    for category, cols in metrics.items():
        if cols:
            for metric_info in cols:
                print(f"  {metric_info['name']}: {category} ({metric_info['granularity']}-level)")

    # Count total metrics to plot
    total_metrics = sum(len(cols) for cols in metrics.values())

    if total_metrics == 0:
        print("No metrics found to plot")
        return

    # Determine subplot layout
    # Priority: loss, accuracy, learning_rate, other
    plot_configs = []

    for metric_info in metrics['loss']:
        plot_configs.append(('loss', metric_info))
    for metric_info in metrics['accuracy']:
        plot_configs.append(('accuracy', metric_info))
    for metric_info in metrics['learning_rate']:
        plot_configs.append(('learning_rate', metric_info))
    for metric_info in metrics['other']:
        plot_configs.append(('other', metric_info))

    # Create subplots
    n_plots = len(plot_configs)

    if n_plots <= 2:
        nrows, ncols = 1, n_plots
    elif n_plots <= 4:
        nrows, ncols = 2, 2
    elif n_plots <= 6:
        nrows, ncols = 2, 3
    else:
        nrows = (n_plots + 2) // 3
        ncols = 3

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    fig.suptitle('Training Process Visualization', fontsize=16, fontweight='bold')

    # Flatten axes for easy iteration
    if n_plots == 1:
        axes = [axes]
    elif nrows > 1 or ncols > 1:
        axes = axes.flatten()
    else:
        axes = [axes]

    # Plot each metric
    for idx, (category, metric_info) in enumerate(plot_configs):
        if idx < len(axes):
            plot_single_metric(data, metric_info, axes[idx])

    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()

    # Save or show
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nFigure saved to {output_path}")
    else:
        plt.show()


def print_summary(data: List[Tuple[str, pd.DataFrame]]):
    """Print summary statistics for each dataset."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)

    for label, df in data:
        print(f"\n{label}:")
        print("-" * 40)

        # Basic info
        print(f"Total steps: {df['step'].max() if 'step' in df.columns else 'N/A'}")
        print(f"Total epochs: {df['epoch'].max() + 1 if 'epoch' in df.columns else 'N/A'}")

        # Metric statistics
        metrics = detect_metric_columns(df)

        for category, metric_infos in metrics.items():
            for metric_info in metric_infos:
                col = metric_info['name']
                granularity = metric_info['granularity']

                if col in df.columns:
                    values = df[col].dropna()
                    if len(values) > 0:
                        print(f"  {col} ({granularity}-level):")
                        print(f"    min: {values.min():.6f}")
                        print(f"    max: {values.max():.6f}")
                        print(f"    mean: {values.mean():.6f}")
                        print(f"    final: {values.iloc[-1]:.6f}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize PyTorch Lightning training process from CSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Visualize single training run
    python visualize_training_process.py training_log.csv

    # Compare multiple runs
    python visualize_training_process.py run1.csv run2.csv run3.csv

    # Save to file
    python visualize_training_process.py training_log.csv -o training_curves.png

    # Custom figure size
    python visualize_training_process.py training_log.csv --figsize 20 12
        """
    )

    parser.add_argument('csv_files', nargs='+', help='CSV file(s) to visualize')
    parser.add_argument('-o', '--output', help='Output path for the figure (PNG, PDF, etc.)')
    parser.add_argument('--figsize', nargs=2, type=int, default=[15, 10],
                       help='Figure size as width height (default: 15 10)')
    parser.add_argument('--no-summary', action='store_true',
                       help='Skip printing summary statistics')
    parser.add_argument('--style', default='seaborn-v0_8-darkgrid',
                       help='Matplotlib style (default: seaborn-v0_8-darkgrid)')

    args = parser.parse_args()

    # Set style
    try:
        plt.style.use(args.style)
    except OSError:
        print(f"Warning: style '{args.style}' not found, using default")
        plt.style.use('seaborn-v0_8-whitegrid')

    # Load data
    data = load_csv_files(args.csv_files)

    if not data:
        print("No valid CSV files loaded")
        sys.exit(1)

    # Print summary
    if not args.no_summary:
        print_summary(data)

    # Create visualization
    plot_training_curves(data, args.output, tuple(args.figsize))


if __name__ == '__main__':
    main()