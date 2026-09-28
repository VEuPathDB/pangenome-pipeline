"""Interval arithmetic.

The pipeline mixes two coordinate conventions, so each function states which one it expects:
  - half-open [start, end), 0-based - PAF/BED-style (minimap2 alignments in step 01)
  - closed [start, end], 1-based    - GFF-style (gene/CDS features)
"""


def merged_len(intervals):
    """Total length covered by half-open [start, end) intervals, counting overlaps once."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs
            cs, ce = s, e
    return total + (ce - cs)


def segments_overlap_bp(segs1, segs2):
    """Overlapping bp between two lists of closed, 1-based [start, end] segments (e.g. two genes' CDS).

    Summed pair by pair, so segments within one list are assumed not to overlap each other.
    """
    total = 0
    for s1, e1 in segs1:
        for s2, e2 in segs2:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi >= lo:
                total += hi - lo + 1
    return total
