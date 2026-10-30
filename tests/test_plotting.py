import pytest
from matplotlib.figure import Figure

from regionreport import annotate_regions, plot_overview, plot_region


def test_feature_cap_is_explicit_and_tables_remain_complete(bed, tmp_path, monkeypatch):
    result = annotate_regions(
        bed([("c", 0, 100, "query")]),
        bed([("c", i, i + 10, f"g{i}") for i in range(10)]),
    )

    def inspect(figure, path, **kwargs):
        assert len(figure.axes[0].patches) == 4  # query + 3 selected genes
        assert "Showing 3 of 10" in figure._supxlabel.get_text()

    monkeypatch.setattr(Figure, "savefig", inspect)
    plot_region(result, "query", tmp_path / "limited.png", max_features=3)
    assert len(result.overlaps) == 10


def test_overview_rejects_overlapping_windows(bed, tmp_path):
    result = annotate_regions(
        bed([("c", 0, 10, "a"), ("c", 5, 20, "b")]),
        bed([]),
    )
    with pytest.raises(ValueError, match="contiguous and non-overlapping"):
        plot_overview(result, tmp_path / "invalid.png")


def test_circular_fragments_share_a_plot_row(bed, tmp_path, monkeypatch):
    result = annotate_regions(
        bed([("c", 90, 10, "query")]),
        bed([("c", 95, 5, "g")]),
        circular=True,
        sequence_lengths={"c": 100},
    )

    def inspect(figure, path, **kwargs):
        axis = figure.axes[0]
        assert len(axis.patches) == 3  # query + 2 pieces of one gene
        assert len(axis.get_yticklabels()) == 2
        assert any(text.get_text() == "sequence coordinate 0" for text in axis.texts)

    monkeypatch.setattr(Figure, "savefig", inspect)
    plot_region(result, "query", tmp_path / "wrapped.svg")
