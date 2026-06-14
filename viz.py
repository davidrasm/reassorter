#!/usr/bin/env python3
"""Draw a simulated ARG as a phylogenetic network with matplotlib.

Consumes the ARGGraph produced by `arg.simulate_arg(..., record_graph=True)`:
a DAG whose nodes are samples, coalescences, and recombination/reassortment
events. Time is drawn on the y-axis (samples at the bottom). Sample and
coalescence nodes are circles; recombination/reassortment nodes are squares
(one child below, two parents above) annotated with their breakpoint (hudson)
or the segments routed to each parent (reassortment).

This is only legible for small ARGs (few samples, modest recombination); large
graphs become an unreadable hairball.

Requires matplotlib (`pip install matplotlib`). The ARGGraph edges already
carry the ancestral segments per lineage, so a future "segments per lineage"
overlay is purely a rendering addition here.
"""


def _seg_indices(intervals):
    """Expand integer-aligned [left, right) intervals into segment indices."""
    idx = []
    for left, right in intervals:
        idx.extend(range(int(round(left)), int(round(right))))
    return idx


def _recomb_label(node):
    meta = node.meta or {}
    if "breakpoint" in meta:
        return f"x={meta['breakpoint']:.3g}"
    if "to_first" in meta:
        first = ",".join(map(str, _seg_indices(meta["to_first"])))
        second = ",".join(map(str, _seg_indices(meta["to_second"])))
        return f"{{{first}}} | {{{second}}}"
    return ""


def _layout(graph):
    """Assign (x, y) to every node: y = time, x = mean of its children's x."""
    children = {}
    for edge in graph.edges:
        children.setdefault(edge.parent, []).append(edge.child)

    samples = sorted(nd.id for nd in graph.nodes if nd.type == "sample")
    x = {sid: float(i) for i, sid in enumerate(samples)}

    # Children are always younger, so processing by increasing time guarantees
    # a node's children already have an x before it is placed.
    for node in sorted(graph.nodes, key=lambda nd: nd.time):
        if node.id in x:
            continue
        kids = children.get(node.id, [])
        x[node.id] = sum(x[c] for c in kids) / len(kids) if kids else 0.0

    y = {nd.id: nd.time for nd in graph.nodes}
    return x, y


def draw_arg(graph, mode="hudson", ax=None, save=None, show=False,
             label_recomb=True):
    """Draw an ARGGraph as a network and return (figure, axes).

    Parameters
    ----------
    graph : arg.ARGGraph
        From simulate_arg(..., record_graph=True).
    mode : str
        Used only for the title ("hudson" or "reassortment").
    ax : matplotlib Axes, optional
        Draw into this axes instead of creating a new figure.
    save : str, optional
        If given, write the figure to this path.
    show : bool
        If True, call plt.show().
    label_recomb : bool
        Annotate recombination nodes with their breakpoint / segment routing.
    """
    if graph is None:
        raise ValueError(
            "graph is None; call simulate_arg(..., record_graph=True)")

    import matplotlib
    if save is not None and not show:
        matplotlib.use("Agg")  # headless: render straight to a file
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    x, y = _layout(graph)
    n_samples = sum(1 for nd in graph.nodes if nd.type == "sample")

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(5.0, 1.3 * n_samples), 6.0))
    else:
        fig = ax.figure

    for edge in graph.edges:
        ax.plot([x[edge.child], x[edge.parent]],
                [y[edge.child], y[edge.parent]],
                color="0.6", lw=1.0, zorder=1)

    styles = {"sample": ("o", "black"),
              "coalescence": ("o", "#1f77b4"),
              "recombination": ("s", "#d62728")}
    for node in graph.nodes:
        marker, color = styles[node.type]
        ax.scatter([x[node.id]], [y[node.id]], marker=marker, color=color,
                   s=45, zorder=3)
        if node.type == "sample":
            ax.annotate(f"n{node.id + 1}", (x[node.id], y[node.id]),
                        textcoords="offset points", xytext=(0, -12),
                        ha="center", fontsize=8)
        elif node.type == "recombination" and label_recomb:
            label = _recomb_label(node)
            if label:
                ax.annotate(label, (x[node.id], y[node.id]),
                            textcoords="offset points", xytext=(7, 0),
                            va="center", fontsize=7, color="#d62728")

    ax.set_ylabel("time (generations before present)")
    ax.set_xticks([])
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.margins(x=0.12, y=0.05)
    ax.set_title("ARG — Hudson recombination" if mode == "hudson"
                 else "ARG — reassortment")

    legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
               markersize=7, label="sample"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f77b4",
               markersize=7, label="coalescence"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#d62728",
               markersize=7, label="recombination"),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=8, frameon=False)

    fig.tight_layout()
    if save is not None:
        fig.savefig(save, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax
