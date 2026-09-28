"""lib/pangenome_common - helpers shared across step scripts."""

from pangenome_common.intervals import merged_len, segments_overlap_bp
from pangenome_common.io import read_strains


def test_merged_len_half_open():
    assert merged_len([]) == 0
    assert merged_len([(0, 10)]) == 10
    assert merged_len([(20, 30), (0, 10), (5, 15)]) == 25  # [0,15) + [20,30)
    assert merged_len([(0, 10), (10, 20)]) == 20  # touching half-open intervals don't double count
    assert merged_len([(0, 100), (10, 20)]) == 100  # nested


def test_segments_overlap_bp_closed():
    assert segments_overlap_bp([(1, 10)], [(10, 20)]) == 1  # 1-based inclusive: position 10 is shared
    assert segments_overlap_bp([(1, 10), (20, 30)], [(5, 25)]) == 12
    assert segments_overlap_bp([(1, 10)], [(11, 20)]) == 0
    assert segments_overlap_bp([], [(1, 5)]) == 0


def test_read_strains_from_file_and_string(tmp_path):
    order = tmp_path / "genome_order.txt"
    order.write_text("# reference first\nREF\n\n  S2  \nS3\n")
    assert read_strains(str(order)) == ["REF", "S2", "S3"]
    assert read_strains("REF, S2,,S3") == ["REF", "S2", "S3"]
