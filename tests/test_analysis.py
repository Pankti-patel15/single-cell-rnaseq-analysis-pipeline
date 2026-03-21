import unittest
from pathlib import Path
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from single_cell_pipeline.analysis import assign_cluster_cell_types, summarize_clusters


class AnnotationTests(unittest.TestCase):
    def test_assign_cluster_cell_types(self):
        obs = pd.DataFrame(
            {
                "leiden": ["0", "0", "1", "1"],
                "T_cells_score": [1.2, 0.9, 0.1, 0.2],
                "B_cells_score": [0.3, 0.2, 1.4, 1.1],
            },
            index=["cell1", "cell2", "cell3", "cell4"],
        )
        annotated = assign_cluster_cell_types(obs)
        self.assertTrue((annotated.loc[["cell1", "cell2"], "cluster_cell_type"] == "T_cells").all())
        self.assertTrue((annotated.loc[["cell3", "cell4"], "cluster_cell_type"] == "B_cells").all())

    def test_summarize_clusters(self):
        obs = pd.DataFrame(
            {
                "leiden": ["0", "0", "1"],
                "cluster_cell_type": ["T_cells", "T_cells", "B_cells"],
            }
        )
        summary = summarize_clusters(obs)
        self.assertEqual(summary["n_cells"].tolist(), [2, 1])


if __name__ == "__main__":
    unittest.main()

