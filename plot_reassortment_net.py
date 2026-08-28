#!/usr/bin/env python3
"""
    Draw a reassortment network as a "subway map": every lineage is drawn as one
    parallel coloured band per segment it carries, so a lineage carrying
    segments 0 and 1 is a blue line and an orange line running side by side.
    Bands start and stop at reassortment events, where the segments they carry
    change.

    Reads extended Newick in a NEXUS wrapper: either a CoalRe summary network
    (dated tips, posterior support) or a simulated ARG from
    `arg.py --mode reassortment --format extended-newick` (tips at t=0, no
    posteriors). The only annotation required on a branch is
    `[&segments={i, j, ...}]`.

        python3 plot_reassortment_net.py sim.tree -o sim-network.png
        python3 plot_reassortment_net.py tswv-summary.tree \
            --segment-order 2 1 0 --segment-names S M L --xmin 1980
"""

import argparse
import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection
from matplotlib.transforms import offset_copy

import baltic as bt


"""
    Lane colours, assigned to lanes in this fixed order and never cycled.

    Lanes sit at a constant offset per segment, so a segment that is absent
    leaves a gap rather than closing it: only consecutive lanes ever touch.
    That makes the adjacent-pair check the right one, and this order passes it
    on a light surface -- worst adjacent CVD deltaE 9.1 (>= 8 target), worst
    normal-vision deltaE 19.6 (>= 15 floor). Aqua, yellow and magenta fall below
    3:1 contrast on white, so the legend is always drawn: identity is carried by
    lane position and the legend, not by colour alone.

    Eight covers the segmented viruses that exist (influenza has 8, TSWV has 3).
    Past eight there is no ninth hue that stays separable, so draw a subset with
    --segment-order rather than adding one.
"""
SEGMENT_COLOURS = ['#2a78d6',  # blue
                   '#eb6834',  # orange
                   '#1baf7a',  # aqua
                   '#eda100',  # yellow
                   '#e87ba4',  # magenta
                   '#008300',  # green
                   '#4a3aa7',  # violet
                   '#e34948']  # red


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Draw a reassortment network with one coloured band per "
                    "segment carried by each lineage.")
    p.add_argument("tree", help="NEXUS file holding an extended-Newick network")
    p.add_argument("-o", "--output", default=None,
                   help="output image; the extension selects the format "
                        "(default: <tree>-network.png)")

    p.add_argument("--segment-order", type=int, nargs="+", default=None,
                   help="segment indices as they appear in the file, listed in "
                        "the order their lanes should be drawn (default: 0, 1, "
                        "... in file order). Also selects a subset: CoalRe "
                        "numbers segments alphabetically, so TSWV's L, M, S "
                        "become biological S, M, L order with '2 1 0'")
    p.add_argument("--segment-names", nargs="+", default=None,
                   help="legend labels for the lanes, in draw order "
                        "(default: seg<i> using the file's own indices)")

    p.add_argument("--time-axis", choices=("auto", "dates", "generations"),
                   default="auto",
                   help="'dates' reads sampling dates out of the tip names; "
                        "'generations' puts the tips at 0 and the root at "
                        "-treeHeight; 'auto' (default) tries dates and falls "
                        "back to generations")
    p.add_argument("--tip-regex", default=r"\_([0-9\-]+)$",
                   help="[--time-axis dates] regex capturing the date in a tip "
                        "name")
    p.add_argument("--date-fmt", default="%Y-%m-%d",
                   help="[--time-axis dates] format of the captured date")
    p.add_argument("--xmin", type=float, default=None,
                   help="left edge of the time axis (default: the root)")

    p.add_argument("--posterior-cutoff", type=float, default=0.01,
                   help="hide reassortments supported below this posterior; "
                        "branches with no posterior annotation count as 1.0 "
                        "(default: 0.01)")
    p.add_argument("--min-branch-length", type=float, default=0.001,
                   help="shortest branch allowed when a summary network has "
                        "negative branch lengths; 0 just removes the "
                        "inversions, and a negative value reports them and "
                        "changes nothing (default: 0.001)")

    p.add_argument("--lane-gap", type=float, default=3.0,
                   help="spacing between neighbouring bands, in points")
    p.add_argument("--lane-width", type=float, default=3.0,
                   help="width of each band, in points")
    p.add_argument("--figsize", type=float, nargs=2, default=None,
                   metavar=("W", "H"),
                   help="figure size in inches, before the side panels are "
                        "added (default: 20 by 0.6 inches per tip)")
    p.add_argument("--dpi", type=int, default=300, help="output DPI")

    p.add_argument("--tip-reassortments", action="store_true",
                   help="label each tip with the number of reassortment events "
                        "on its path from the root")
    p.add_argument("--no-tip-reassortment-bars", dest="tip_reassortment_bars",
                   action="store_false",
                   help="drop the bar panel on the right showing those same "
                        "counts, one bar per tip")
    p.add_argument("--segment-boxes", action="store_true",
                   help="draw rectangles in the right margin showing which "
                        "segments moved at each event")
    p.add_argument("--posterior-labels", action="store_true",
                   help="write each event's posterior support in the right "
                        "margin")
    p.add_argument("--no-ltt", dest="ltt", action="store_false",
                   help="drop the panel below counting the lineages carrying "
                        "each segment through time")
    p.add_argument("--no-flush-tips", dest="flush_tips", action="store_false",
                   help="let each segment lane of an external branch end "
                        "staggered by its lane offset instead of squared off "
                        "on the sampling time")
    return p.parse_args(argv)


def load_network(args):
    """Load the NEXUS network and place it in time. Returns (tree, dated)."""
    if args.time_axis in ("auto", "dates"):
        try:
            return bt.loadNexus(args.tree, tip_regex=args.tip_regex,
                                date_fmt=args.date_fmt), True
        except (AssertionError, ValueError):
            if args.time_axis == "dates":
                raise SystemExit(
                    "no sampling dates found in the tip names with "
                    f"--tip-regex {args.tip_regex!r} and --date-fmt "
                    f"{args.date_fmt!r}; use --time-axis generations")

    ll = bt.loadNexus(args.tree, absoluteTime=False)
    ll.setAbsoluteTime(0.0)  # tips at 0, root at -treeHeight
    return ll, False


def main(argv=None):
    args = parse_args(argv)
    output = args.output or (os.path.splitext(args.tree)[0] + "-network.png")

    ll, dated = load_network(args)

    """
        A reassortment is one instant reached by two different paths -- the
        reticulation branch and the node it lands on -- and everything below
        relies on the two agreeing exactly, so the jump edges stay vertical and
        the lineage counts balance. Newick carries branch lengths as decimal
        text, so the two paths accumulate rounding differently and can disagree
        in the last digit even on an exactly simulated network. Snap them
        together once, up front, rather than depending on the file's precision.
    """
    for k in ll.Objects:
        if isinstance(k, bt.reticulation):
            t = max(k.absoluteTime, k.target.absoluteTime)
            k.absoluteTime = k.target.absoluteTime = t

    """
        Which segments exist is read off the network itself rather than assumed:
        a branch's [&segments={...}] annotation lists the segments it carries,
        and baltic parses those into floats.
    """
    present = sorted({int(s) for k in ll.Objects
                      for s in k.traits.get('segments', [])})
    if not present:
        raise SystemExit(f"no [&segments={{...}}] annotations in {args.tree}; "
                         "this is not a segment-annotated network")

    order = args.segment_order if args.segment_order is not None else present
    unknown = [s for s in order if s not in present]
    if unknown:
        raise SystemExit(f"--segment-order names segment(s) {unknown} that do "
                         f"not appear in {args.tree} (it carries {present})")
    if len(order) > len(SEGMENT_COLOURS):
        raise SystemExit(f"{len(order)} segments exceeds the {len(SEGMENT_COLOURS)} "
                         "lane colours that stay separable; draw a subset with "
                         "--segment-order")

    segment_map = {seg: lane for lane, seg in enumerate(order)}  # file index -> lane
    n_lanes = len(order)
    segment_names = args.segment_names or [f"seg{seg}" for seg in order]
    if len(segment_names) != n_lanes:
        raise SystemExit(f"--segment-names has {len(segment_names)} labels for "
                         f"{n_lanes} lanes")
    segment_cmap = {lane: SEGMENT_COLOURS[lane] for lane in range(n_lanes)}

    lane_gap = args.lane_gap
    lane_width = args.lane_width

    def carried(k):
        """Lanes carried by branch k (empty for the root stub)."""
        return {segment_map[int(s)] for s in k.traits.get('segments', [])
                if int(s) in segment_map}

    """
        Simulated networks carry no posterior support, so a missing posterior
        reads as full support: every event is drawn, at full opacity, and the
        support-related decorations switch themselves off.
    """
    has_posterior = any('posterior' in k.traits for k in ll.Objects)
    posterior = lambda k: k.traits.get('posterior', 1.0)
    posteriorCutoff = lambda k: posterior(k) >= args.posterior_cutoff

    alpha_floor = 0.1  # opacity of a posterior 0 event

    def edge_alpha(k):
        """Opacity for a reassortment event, scaled by its posterior support."""
        return alpha_floor + (1 - alpha_floor) * min(max(posterior(k), 0.0), 1.0)

    show_posterior_labels = args.posterior_labels and has_posterior

    def inverted(ll):
        """Branches that start before the parent they descend from."""
        return [k for k in ll.Objects
                if k.parent is not None and k.parent.absoluteTime is not None
                and k.absoluteTime < k.parent.absoluteTime - 1e-9]

    units = 'years' if dated else 'generations'
    negative = [k for k in ll.Objects
                if k.length is not None and k.length < 0]  # summary networks can carry these
    if negative:
        print('%d branch(es) of %d have negative length (worst %.4f), putting '
              '%d child branch(es) before their parent' % (
                  len(negative), len(ll.Objects),
                  min(k.length for k in negative), len(inverted(ll))))

    if negative and args.min_branch_length >= 0:
        """
            Two passes run to a fixed point.

            The sweep pushes any branch starting before its parent forward to
            parent + min_branch_length. It works on absoluteTime rather than
            clamping length and re-deriving the times, and that choice matters:
            on a summary network, clamping lengths and recomputing heights
            rigidly shifts whole subtrees and drags sampled tips off their
            collection dates. Sweeping the times moves only what has to move, so
            tips stay where the data put them.

            The sweep alone is not enough. A reticulation and the node it lands
            on are one event reached by two different paths, so moving one path
            pulls them apart -- which would tilt the jump edges off vertical and
            break the lineage counts below, which rely on an event ending one
            branch and starting another at the same instant. Each pair is pulled
            back onto the later of its two times, and the sweep runs again in
            case that push created a fresh inversion. Times only ever increase,
            so this settles.
        """
        for _ in range(100):
            moved = False

            stack = [ll.root]
            while stack:  # no branch may begin before its parent
                k = stack.pop()
                if k.parent is not None and k.parent.absoluteTime is not None:
                    floor = k.parent.absoluteTime + args.min_branch_length
                    if k.absoluteTime < floor:
                        k.absoluteTime = floor
                        moved = True
                if k.is_node():
                    stack.extend(k.children)

            for k in ll.Objects:  # a reassortment happens at one instant, on both paths reaching it
                if not isinstance(k, bt.reticulation):
                    continue
                t = max(k.absoluteTime, k.target.absoluteTime)
                if k.absoluteTime != t or k.target.absoluteTime != t:
                    moved = True
                k.absoluteTime = k.target.absoluteTime = t

            if not moved:
                break
        else:
            raise RuntimeError('branch length correction did not settle')

        for k in ll.Objects:  # bring lengths back in line with the times now being drawn
            if k.parent is not None and k.parent.absoluteTime is not None:
                k.length = k.absoluteTime - k.parent.absoluteTime
        ll.traverse_tree()  # re-derive heights from the corrected lengths (leaves the y layout alone)

        assert not inverted(ll), 'branches still run backwards after correction'
        assert max((abs(k.absoluteTime - k.target.absoluteTime)
                    for k in ll.Objects if isinstance(k, bt.reticulation)),
                   default=0.0) < 1e-9, 'reassortments no longer land at one instant'
        print('  corrected to a minimum branch length of %g %s'
              % (args.min_branch_length, units))

    tips = [k for k in ll.Objects if isinstance(k, bt.leaf)]  # sampled tips only -- reticulations also carry branchType 'leaf'

    """
        The network keeps its full size whatever else is switched on; the extra
        inches are added to the figure rather than taken out of the tree. Each
        side panel is tied to the network by a shared axis rather than by
        matching coordinates up by hand -- x for the lineage panel below, so it
        can be read straight down from any point in time, and y for the bar
        panel to the right, so every bar sits on the row of the tip it counts.
    """
    width, height = args.figsize or (20.0, max(8.0, 0.6 * len(tips)))
    fig = plt.figure(figsize=(width + (3 if args.tip_reassortment_bars else 0),
                              height + (height / 4 if args.ltt else 0)),
                     facecolor='w')
    gs = gridspec.GridSpec(2 if args.ltt else 1,
                           2 if args.tip_reassortment_bars else 1,
                           height_ratios=[4, 1] if args.ltt else [1],
                           width_ratios=[20, 2.5] if args.tip_reassortment_bars else [1],
                           hspace=0.06, wspace=0.02)
    ax = plt.subplot(gs[0, 0], facecolor='w')  # the network
    ax_ltt = plt.subplot(gs[1, 0], facecolor='w', sharex=ax) if args.ltt else None  # lineages through time, same time axis
    ax_bar = plt.subplot(gs[0, 1], facecolor='w', sharey=ax) if args.tip_reassortment_bars else None  # reassortment counts, same tips

    mostRecentTip = max(ll.getParameter('absoluteTime', use_trait=False))

    """
        Count the reassortment events on each branch's ancestry, root to tips.

        baltic hangs a .contribution on the node a reticulation lands on, so
        those nodes are the events themselves; walking down from the root with a
        running total leaves every branch holding the number of events above it.
        Events are counted subject to the same posterior cutoff used to draw
        them, so the numbers agree with what is on screen.

        This follows the drawn topology. In the network a lineage also traces
        back through every reticulation it received, so a tip's full ancestry is
        a DAG with more than one path to the root -- what is counted here is the
        single path baltic draws, the one the eye can actually follow up the
        figure.
    """
    stack = [(ll.root, 0)]
    while stack:
        w_branch, n = stack.pop()
        if hasattr(w_branch, 'contribution') and posteriorCutoff(w_branch.contribution):
            n += 1  # a reassortment event lands here
        w_branch.traits['n_reassort'] = n
        if w_branch.is_node():
            stack.extend((c, n) for c in w_branch.children)  # carry the running total down

    """
        The right margin only holds whatever is switched on, so it collapses
        rather than leaving an empty gutter, and the axis stays tight in every
        combination. The dashed leader line exists purely to carry the eye from
        an event out to its row in that margin, so with the margin empty it has
        nothing to point at and is dropped along with it.

        The limits are settled here, before anything is drawn, because
        converting the lane offsets from points into time units below needs a
        final x axis. Everything in the margin is sized as a fraction of the
        time span rather than in fixed units, so it looks the same whether the
        axis runs over decades or thousands of generations.
    """
    show_leader_lines = args.segment_boxes or show_posterior_labels

    x_left = args.xmin if args.xmin is not None else min(
        ll.getParameter('absoluteTime', use_trait=False))
    span = mostRecentTip - x_left
    w = 0.0079 * span  # width of an indicator box
    h = 2              # height of a box, in tip rows
    step = 0.0214 * span   # spacing unit out in the margin
    base = mostRecentTip + 0.5 * step  # where the indicator boxes begin

    box_span = n_lanes * w if args.segment_boxes else 0.0
    label_x = base + box_span + step  # where the posterior labels sit
    if show_posterior_labels:
        xmax = label_x + 3 * step  # room for the label text
    elif args.segment_boxes:
        xmax = base + box_span + step  # just the boxes
    else:
        xmax = base + step  # nothing out here, so stop just past the tips

    ax.set_xlim(x_left - 0.02 * span, xmax)

    """
        One transform per segment, offsetting that segment's whole lane by a
        constant vector measured in points. Because the shift is uniform across
        the lane, horizontal branches separate vertically, vertical node bars
        separate horizontally, and every corner still meets exactly. Working in
        points (rather than data units) keeps the spacing between bands constant
        regardless of figure size or axis limits.
    """
    lane_offset = {}
    lane_transform = {}
    for s in range(n_lanes):
        lane_offset[s] = ((n_lanes - 1) / 2 - s) * lane_gap  # first segment rides on top
        lane_transform[s] = offset_copy(ax.transData, fig=fig, x=lane_offset[s],
                                        y=lane_offset[s], units='points')

    def units_per_point():
        """One typographic point, expressed in time units on the current axis."""
        inv = ax.transData.inverted()
        return inv.transform((fig.dpi / 72.0, 0))[0] - inv.transform((0, 0))[0]  # dpi cancels out, so this holds at any savefig resolution

    tip_pullback = units_per_point() if args.flush_tips else 0.0  # for squaring off the tips, computed once

    for s in range(n_lanes):  # draw one lane per segment
        lines = []
        for k in ll.Objects:
            if isinstance(k, bt.reticulation):
                continue  # reassorting branches are drawn dashed further down

            """
                The branch above a node and the bar below it are decided
                separately, because a node can need one without the other. A
                segment is dropped from a lineage once it has reached its MRCA,
                so the branch above that last coalescence does not carry it (and
                the root branch carries nothing at all) -- but the coalescence
                itself still has to be drawn, or the segment's two final lineages
                run to the same time and stop, appearing never to meet.
            """
            here = s in carried(k)
            below = [c.y for c in k.children if s in carried(c)] if k.is_node() else []
            if not here and len(below) < 2:
                continue  # this segment neither passes through here nor coalesces here

            x = k.absoluteTime
            y = k.y

            if here:
                xp = k.parent.absoluteTime if k.parent else x

                end = x
                if args.flush_tips and isinstance(k, bt.leaf):
                    """
                        Every lane is shifted sideways by its own offset, which
                        leaves the lanes of a tip branch ending a few points
                        apart. Pull the endpoint back by exactly that offset so
                        the lane's own shift carries it to the sampling time and
                        the tip squares off. Only the far end moves; the near end
                        keeps its offset and stays joined to the vertical bar
                        above it.
                    """
                    end = x - lane_offset[s] * tip_pullback

                lines.append(((xp, y), (end, y)))  # horizontal branch

            if below:  # vertical bar joining this node to its descendants
                """
                    Only span children that actually carry this segment -- at a
                    reassortment the other child inherits a disjoint set, and a
                    bar drawn to children[0]..children[-1] (what baltic's
                    plotTree does) would run a band into a lineage that never
                    carries it. The node's own y joins in only when the segment
                    continues above, which keeps the bar attached when every
                    carrying child lies to one side; at an MRCA there is nothing
                    above to attach to and the bar just spans the two children.
                """
                ys = below + ([y] if here else [])
                if len(ys) > 1:
                    lines.append(((x, min(ys)), (x, max(ys))))

        ax.add_collection(LineCollection(lines, lw=lane_width, color=segment_cmap[s],
                                         transform=lane_transform[s],
                                         capstyle='projecting', zorder=10 + s))

    for k in ll.Objects:
        if isinstance(k, bt.reticulation) and posteriorCutoff(k):  # a reassorting branch
            size = 120  # size of the triangle marking the event
            xp = k.parent.absoluteTime
            x = k.absoluteTime
            y = k.y

            segs = carried(k)  # travelling segment(s), in lane order
            a = edge_alpha(k)  # fade the whole event -- edge, marker and indicator -- by its support
            zo = posterior(k)  # and let better supported events draw over weaker ones

            for s in sorted(segs):  # a reassortment can carry more than one segment
                c = segment_cmap[s]
                tr = lane_transform[s]  # draw in that segment's lane so the jump leaves and lands on the right band

                if show_leader_lines:
                    ax.plot([x, base], [y, y], ls='--', color=c, lw=1, alpha=a,
                            zorder=2 + zo, transform=tr)  # from the reassorting lineage out to the margin

                ax.plot([xp, x, x], [y, y, k.target.y], color=c, lw=lane_width,
                        ls='--', alpha=a, zorder=5 + zo, transform=tr)  # leave the donor, head to the recipient

                marker = '^' if k.target.y > k.y else 'v'  # point towards the recipient
                ax.scatter(x, k.target.y, s=size, marker=marker, facecolor=c,
                           edgecolor='w', lw=1, alpha=a, zorder=300 + zo, transform=tr)

            if show_posterior_labels:
                ax.text(label_x, y, '%.2f' % posterior(k), size=14, color='k',
                        alpha=a, ha='left', va='center')

            if args.segment_boxes:
                for s in range(n_lanes):
                    fc = segment_cmap[s] if s in segs else 'w'  # coloured if this segment is reassorting
                    ax.add_patch(plt.Rectangle((base + s * w, y - h / 2), w, h,
                                               facecolor=fc, edgecolor='k', lw=1,
                                               alpha=a, clip_on=False))

    if args.tip_reassortments:  # number each tip with the reassortments on its path from the root
        tip_offset = offset_copy(ax.transData, fig=fig, x=8, y=0, units='points')  # sit just clear of the tip
        for k in tips:
            ax.text(k.absoluteTime, k.y, '%d' % k.traits['n_reassort'],
                    transform=tip_offset, size=11, color='0.25', ha='left',
                    va='center', zorder=400)

    if args.tip_reassortment_bars:  # the same counts as bars, each on the row of its tip
        """
            Sharing y with the network is what lines the bars up: the panels
            hold the same y coordinates and the same height, so a tip and its bar
            land on one row without any positions being matched up by hand. Tips
            sit 1 to 1.5 apart in y (reticulations take up slots between them),
            so a bar height of 0.8 clears its neighbours at the tightest spacing.
        """
        counts = [k.traits['n_reassort'] for k in tips]
        ax_bar.barh([k.y for k in tips], counts, height=0.8, color='0.35', zorder=3)

        [ax_bar.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
        ax_bar.set_xlim(0, max(counts) + 0.5)
        ax_bar.set_xticks(range(0, max(counts) + 1))
        ax_bar.tick_params(labelsize=20)
        ax_bar.set_xlabel('reassortments', size=20)
        ax_bar.xaxis.grid(True, ls=':', lw=1, color='0.8')  # recessive gridlines, behind the bars
        ax_bar.set_axisbelow(True)

    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    ax.xaxis.grid(True, ls='--', lw=1, color='grey')  # time lines, wherever the ticks land
    ax.set_axisbelow(True)
    ax.set_yticks([])
    ax.set_yticklabels([])

    ax.tick_params(labelsize=20)
    ax.set_ylim(-1, ll.ySpan + 1)  # the x limits were fixed before drawing, so the lane offsets could be converted into time units

    time_label = 'date' if dated else 'generations before present'

    if args.ltt:  # count the lineages carrying each segment through time
        """
            Counted by interval rather than by event type: every branch
            contributes +1 at the time it starts and -1 at the time it ends, once
            for each segment it carries, and the running sum is the number of
            lineages carrying that segment.

            Doing it this way gets reassortments right for free. An event ends
            one parent branch and starts the child at the same instant, so a
            segment's count steps only at coalescences and samplings and never at
            a reassortment -- which is correct, since a reassortment repackages
            which lineage carries a segment without changing how many carry it.

            The whole network is counted whatever the posterior cutoff is set to.
            Dropping a reticulation would end a lineage with nothing replacing
            it, leaving its segments to reappear out of nowhere further down.
        """
        ltt_max = 0
        for s in range(n_lanes):
            delta = {}
            for k in ll.Objects:
                if s not in carried(k):
                    continue  # this lineage does not carry this segment
                start = round(k.parent.absoluteTime if k.parent else k.absoluteTime, 6)  # rounded so both ends of a reassortment land on one key
                end = round(k.absoluteTime, 6)
                delta[start] = delta.get(start, 0) + 1  # a lineage carrying this segment begins
                delta[end] = delta.get(end, 0) - 1      # ...and ends

            times = sorted(delta)
            counts = []
            n = 0
            for t in times:  # running sum over the sorted event times gives the step function
                n += delta[t]
                counts.append(n)
            ltt_max = max(ltt_max, max(counts))

            offset = ((n_lanes - 1) / 2 - s) * lane_gap  # same lane offsets as the network above, so segments carried by equal numbers of lineages stay legible instead of overplotting
            ax_ltt.step(times, counts, where='post', color=segment_cmap[s],
                        lw=lane_width,
                        transform=offset_copy(ax_ltt.transData, fig=fig,
                                              x=offset, y=offset, units='points'))

        [ax_ltt.spines[loc].set_visible(False) for loc in ['top', 'right']]
        ax_ltt.xaxis.grid(True, ls='--', lw=1, color='grey')  # same time lines as the network
        ax_ltt.set_axisbelow(True)
        ax_ltt.set_ylim(0, ltt_max * 1.12)
        ax_ltt.set_ylabel('lineages', size=20)
        ax_ltt.set_xlabel(time_label, size=20)
        ax_ltt.tick_params(labelsize=20)
        ax.tick_params(labelbottom=False)  # the time axis is labelled once, on the panel below
    else:
        ax.set_xlabel(time_label, size=20)

    handles = [plt.Line2D([], [], color=segment_cmap[s], lw=5, label=segment_names[s])
               for s in range(n_lanes)]  # legend keyed to lane order, top to bottom
    ax.legend(handles=handles, fontsize=20, frameon=False, loc='lower left')  # bottom left is empty; upper left collides with the deep root branches
    if has_posterior:
        ax.text(0.012, 0.135, 'reassortment opacity scales with posterior support',
                transform=ax.transAxes, size=16, color='0.35', va='bottom')  # opacity is an encoding, so say what it means

    plt.savefig(output, dpi=args.dpi, bbox_inches='tight')
    print('wrote %s' % output)


if __name__ == '__main__':
    main()
