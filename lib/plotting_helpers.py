"""Functions to help plotting data from the experiment collection."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from collections import defaultdict
from matplotlib import patches

def plot_single_contig(exp, contig_name, pdf=None):
    """Using data from an experiment collection `exp`, plot a contig."""
    # Get the set of clusters (protein-coding genes) on this contig
    contig_gene_df = exp.contig_df(contig_name)

    # Make the position numeric
    contig_gene_df["start"] = contig_gene_df["start"].apply(int)
    contig_gene_df["end"] = contig_gene_df["end"].apply(int)

    # Calculate the middle of each gene
    contig_gene_df["middle"] = contig_gene_df.loc[
        :,
        ["start", "end"]
    ].mean(axis=1)

    # Sort by length
    contig_gene_df.sort_values(by="start", inplace=True)
    contig_gene_df.reset_index(drop=True, inplace=True)

    gene_colors = dict(zip(contig_gene_df["cluster"].tolist(
    ), sns.color_palette("colorblind", contig_gene_df.shape[0])))

    # Set up the plot
    _, ax = plt.subplots(figsize=[contig_gene_df.shape[0] / 5, 3])

    # Plot the contig
    contig_len = int(contig_name.split("_")[3])

    ax.plot([0, contig_len], [0, 0], c="black", lw=1)

    # Plot each gene
    for ix, r in contig_gene_df.iterrows():
        if r["strand"] == "+":
            height = 0.15
        else:
            height = -0.15

        rect = patches.Rectangle(
            [r["start"], 0],
            r["end"] - r["start"],
            height,
            linewidth=0,
            facecolor=gene_colors[r["cluster"]]
        )
        ax.add_patch(rect)

        # Write the gene name
        text_x_pos = contig_len * ix / contig_gene_df.shape[0]
        plt.text(
            text_x_pos,
            0.65,
            "{} - {}".format(
                r["cluster"], 
                r["product"].rstrip("\n")
            ),
            rotation="vertical",
            verticalalignment="bottom",
            horizontalalignment="center",
        )
        # Connect the gene name to the gene
        ax.plot(
            [r["middle"], r["middle"], text_x_pos, text_x_pos],
            [height, 0.35, 0.5, 0.6],
            color=gene_colors[r["cluster"]],
            lw=1,
        )

    plt.xlabel("Contig position (bp)")
    plt.ylim(-0.5, 1)
    plt.yticks([])

    plt.tight_layout()
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")
    plt.show()


def plot_contig_structure(exp, central_contig_name, min_prop=0.5):
    """Plot the structure of this contig, as compared to any other similar contigs."""

    # Get the set of protein-coding genes on this contig
    clusters = exp.contig_df(central_contig_name)["cluster"].unique()

    # For each cluster, get the set of contigs it was found in (>=90% sequence identity)
    clusters_contigs = {
        cluster_name: exp.contigs_with_gene(cluster_name)
        for cluster_name in clusters
    }

    # Count up the proportion of the clusters found by each contig
    contig_cluster_counts = defaultdict(int)
    for cluster_id, contig_list in clusters_contigs.items():
        for contig_id in contig_list:
            contig_cluster_counts[contig_id] += 1
    contig_cluster_counts = pd.Series(
        contig_cluster_counts).sort_values(ascending=False)

    # Take the proportion
    contig_cluster_prop = contig_cluster_counts / clusters.shape[0]

    # Calculate the proportion of genes that are found
    contigs_to_plot = contig_cluster_prop.index.values[contig_cluster_prop >= min_prop]

    # Get the genes for each contig
    contig_genes = {
        contig_id: exp.contig_df(contig_id)
        for contig_id in contigs_to_plot
    }

    # Get the sample name for each contig
    contig_samples = {
        contig_id: "_".join(contig_gene_df["ID"].values[0].split("_")[:3])
        for contig_id, contig_gene_df in contig_genes.items()
    }

    # Calculate the relative plotting location for the genes in each contig
    central_contig_cluster_start = contig_genes[central_contig_name].set_index("cluster")[
        "start"].apply(int)
    for contig_name, contig_gene_df in contig_genes.items():
        # Convert the coordinates to numbers
        for k in ["start", "end"]:
            contig_gene_df[k] = contig_gene_df[k].apply(int)
            contig_gene_df[k + "_adj"] = contig_gene_df[k]

        # Don't adjust the location of the central contig
        if contig_name == central_contig_name:
            continue

        # Calculate the offset for each individual protein
        offset_df = pd.DataFrame([
            {
                "cluster_id": r["cluster"],
                "position": r["start"],
                "central_contig_position": central_contig_cluster_start[r["cluster"]],
                "offset": r["start"] - central_contig_cluster_start[r["cluster"]],
            }
            for _, r in contig_gene_df.iterrows()
            if r["cluster"] in clusters
        ])
        offset_df.sort_values(by="position", inplace=True)

        # Calculate the marginal offset, comparing one gene to the next immediately adjacent
        marginal_offset = pd.Series(
            offset_df["offset"].values[:-1] - offset_df["offset"].values[1:])

        # If the median marginal offset == 0, the contig is in the correct orientation
        # Otherwise, reverse the orientation of the contig
        if np.abs(marginal_offset.median()) >= 1.:
            contig_gene_df["start"], contig_gene_df["end"] = -1 * \
                contig_gene_df["end"], -1 * contig_gene_df["start"]

            contig_gene_df["strand"] = contig_gene_df["strand"].apply(
                lambda x: "-" if x == "+" else "+")

        # Now calculate the average offset
        avg_offset = np.mean([
            r["start"] - central_contig_cluster_start[r["cluster"]]
            for _, r in contig_gene_df.iterrows()
            if r["cluster"] in clusters
        ])

        # Now calculate the adjusted start and end positions
        contig_gene_df["start_adj"] = contig_gene_df["start"] - avg_offset
        contig_gene_df["end_adj"] = contig_gene_df["end"] - avg_offset

    # Set a color for every gene
    genes = contig_genes[central_contig_name].sort_values(by="start")[
        "cluster"].tolist()
    for gene_df in contig_genes.values():
        genes.extend(
            gene_df.loc[~gene_df["cluster"].isin(genes)].sort_values(by="start")[
                "cluster"].tolist()
        )
    gene_colors = dict(zip(genes, sns.color_palette("colorblind", len(genes))))

    # Set the y coordinates for each contig
    contig_y_coord = {
        contig_name: -ix for ix, contig_name in enumerate(contigs_to_plot)
    }

    # Set up the plot
    fig, ax = plt.subplots(figsize=[10, len(contigs_to_plot) * 0.75])

    # Plot semi-transparent boxes on a per-gene basis
    for cluster_name, cluster_coordinates in pd.concat(list(contig_genes.values())).groupby("cluster"):
        if cluster_name not in clusters:
            continue

        # Assign the y coordinates
        cluster_coordinates["plot_y"] = cluster_coordinates["seqname"].apply(
            contig_y_coord.get)

        # Sort by Y
        cluster_coordinates.sort_values(by="plot_y", inplace=True)

        # Make the polygon
        polygon = patches.Polygon(
            [
                [r["start_adj"], r["plot_y"]]
                for _, r in cluster_coordinates.iterrows()
            ] + [
                [r["end_adj"], r["plot_y"]]
                for _, r in cluster_coordinates.iterrows()
            ][::-1],
            facecolor=gene_colors[cluster_name],
            alpha=0.2,
            closed=True,
            lw=0
        )
        ax.add_patch(polygon)

    # Now plot every contig
    for plot_ix, contig_name in enumerate(contigs_to_plot):
        contig_gene_df = contig_genes[contig_name]

        # Plot the contig itself
        abs_min = contig_gene_df.loc[:, ["start_adj", "end_adj"]].min().min()
        abs_max = contig_gene_df.loc[:, ["start_adj", "end_adj"]].max().max()
        ax.plot([abs_min, abs_max], [-plot_ix, -plot_ix], c="black", lw=1)

        # Plot each gene
        for _, r in contig_gene_df.iterrows():
            if r["strand"] == "+":
                height = 0.25
            else:
                height = -0.25

            rect = patches.Rectangle(
                [r["start_adj"], -plot_ix],
                r["end_adj"] - r["start_adj"],
                height,
                linewidth=0,
                facecolor=gene_colors[r["cluster"]]
            )
            ax.add_patch(rect)

    # Limit the size of the plot to the contig of interest
    central_contig_length = contig_genes[central_contig_name].loc[:, [
        "start", "end"]].max().max()
    plt.xlim(0 - (central_contig_length * 0.1), central_contig_length * 1.1)
    plt.yticks(
        [-1 * x for x in range(len(contigs_to_plot))],
        contigs_to_plot
    )
    plt.xlabel("Contig position (bp)")

    plt.show()
