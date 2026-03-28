import os
import argparse
import numpy as np


class GraphGen:
    def __init__(
            self,
            symmetrized: bool = False,
            weighted: bool = True,
            dense: bool = False,
            seed: int = 42,
    ):
        self.symmetrized = symmetrized
        self.weighted = weighted
        self.dense = dense
        self.rng = np.random.default_rng(seed)

    def _generate(self, n: int) -> np.ndarray:
        # target edge count
        edge_count = (1.0 + self.rng.random()) * n
        if self.dense:
            edge_count *= n

        if self.symmetrized:
            total_possible = n * (n - 1) // 2
        else:
            total_possible = n * (n - 1)

        edge_count = min(edge_count, total_possible)
        p = edge_count / total_possible if total_possible > 0 else 0.0

        if self.symmetrized:
            # Build only upper triangle, then mirror
            mask = self.rng.random((n, n)) < p
            mask = np.triu(mask, k=1)

            if self.weighted:
                weights = self.rng.integers(1, n + 1, size=(n, n))
            else:
                weights = np.ones((n, n), dtype=np.int64)

            adj = np.where(mask, weights, 0)
            adj = adj + adj.T
        else:
            # Directed graph without self-loops
            mask = self.rng.random((n, n)) < p
            np.fill_diagonal(mask, False)

            if self.weighted:
                weights = self.rng.integers(1, n + 1, size=(n, n))
            else:
                weights = np.ones((n, n), dtype=np.int64)

            adj = np.where(mask, weights, 0)

        return adj

    def write(self, n: int, output_folder: str):
        adj = self._generate(n)
        os.makedirs(output_folder, exist_ok=True)
        fl_name = os.path.join(output_folder, f"{n}.txt")

        rows, cols = np.nonzero(adj)
        with open(fl_name, "w") as file:
            for u, v in zip(rows, cols):
                file.write(f"+ {u} {v} {adj[u, v]}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Generate a graph for SSSP benchmarks."
    )

    parser.add_argument(
        "-n",
        type=int,
        default=100,
        help="Number of vertices"
    )

    parser.add_argument(
        "-d", "--dense",
        action="store_true",
        help="Generate a dense graph. Otherwise generate a sparse graph with O(n) edges."
    )

    parser.add_argument(
        "-s", "--symmetrize",
        action="store_true",
        help="Generate an undirected graph by symmetrizing edges."
    )

    parser.add_argument(
        "-w", "--weighted",
        action="store_true",
        help="Generate weighted edges. Otherwise use unweighted edges."
    )

    parser.add_argument(
        "-o", "--output",
        default="data",
        help="Output folder for the graph"
    )

    args = parser.parse_args()

    graph_gen = GraphGen(
        symmetrized=args.symmetrize,
        weighted=args.weighted,
        dense=args.dense
    )
    graph_gen.write(n=args.n, output_folder=args.output)


if __name__ == '__main__':
    main()