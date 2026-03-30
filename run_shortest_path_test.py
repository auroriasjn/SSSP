#!/usr/bin/env python3
import subprocess
import argparse
import shutil
import sys
import re
import os

from datetime import datetime
from dataclasses import dataclass
from typing import Optional

TIME_PATTERN = re.compile(r"Time:\s*([0-9]*\.?[0-9]+)\s*sec", re.IGNORECASE)
MISMATCH_PATTERN = re.compile(r"Total mismatches:\s*(\d+)", re.IGNORECASE)


@dataclass(frozen=True)
class SSSPResults:
    mismatches: int
    seconds: float


# https://github.com/ucrparlay/Parallel-SSSP
def compile_cpp(rebuild=False):
    """
    Compile using Makefile.
    """
    if not shutil.which("make"):
        print("Error: 'make' not found.")
        sys.exit(1)

    try:
        if rebuild:
            subprocess.run(["make", "clean"], check=True)

        subprocess.run(["make"], check=True)
        print("Build successful.\n")

    except subprocess.CalledProcessError:
        print("Build failed.")
        sys.exit(1)


def parse_cpp_output(output: str) -> Optional[SSSPResults]:
    # Extract all solver times (ignore reference times automatically)
    times = [float(t) for t in TIME_PATTERN.findall(output)]

    total_time = float(sum(times)) if times else 0.0

    mismatches = 0
    mismatch_match = MISMATCH_PATTERN.search(output)
    if mismatch_match:
        mismatches = int(mismatch_match.group(1))

    return SSSPResults(mismatches=mismatches, seconds=total_time)


def get_graph_statistics(input_file: str):
    with open(input_file, 'r') as f:
        # The format of the graphs is assumed to be '+ u v w'. Ignore all deletions because
        # this is not batch dynamic.
        edges = [line.strip() for line in f if line.startswith('+')]
        num_edges = len(edges)
        vertices = set()
        for edge in edges:
            _, u, v, _ = edge.split()
            vertices.add(u)
            vertices.add(v)
        num_vertices = len(vertices)
    return num_vertices, num_edges


def run_cpp_executable(
        N: int,
        input_file: str = 'data/usa_batch',
        output_file: str = 'output/output.txt',
        executable: str = "build/sssp_solver",
        algorithm: str = "dijkstra",
        unweighted: bool = False
) -> Optional[SSSPResults]:
    """
    Pack the data in binary (N as int64) and run the C++ executable.
    """
    args = [executable,
            "-i", input_file,
            "-o", output_file,
            "-a", algorithm,
            "-v",
            "-n", str(N)]
    if not unweighted:
        args.append("-w")

    subprocess.run(args, stdout=subprocess.DEVNULL)

    # Parsing the final output
    with open(output_file) as f:
        output = f.read()
        results = parse_cpp_output(output)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run SSSP benchmark tests."
    )

    parser.add_argument(
        "-i", "--input",
        default='data/usa_batch',
        help="Input file for SSSP executable"
    )

    parser.add_argument(
        "-o", "--output",
        default='output',
        help="Output folder for SSSP executable"
    )

    parser.add_argument(
        "-e", "--executable",
        default="build/sssp_solver",
        help="Path to executable (default: ./sssp)"
    )

    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Force clean rebuild before running"
    )

    parser.add_argument(
        "-uw", "--unweighted",
        action="store_true",
        help="Force unweighted graph"
    )

    parser.add_argument(
        "--flush_output",
        action="store_true",
        help="Force flush output folder"
    )

    args = parser.parse_args()
    curr_time = datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
    if args.flush_output and os.path.exists(args.output):
        print("Flushing output directory...")
        shutil.rmtree(args.output)

    if not os.path.exists(args.output):
        print("Making output directory...")
        os.makedirs(args.output)

    if args.unweighted:
        print("Making graph unweighted...")

    # Each run saved in a separate directory
    os.makedirs(f"{args.output}/{curr_time}")
    compile_cpp(rebuild=args.rebuild)

    test_sizes = [10, 100, 1000]
    algorithms = ["dijkstra", "bf", "parallel-dijkstra", "bd", "parallel-bd", "rho-stepping"]
    all_passed = True

    # Iterate over all options
    print(f"Analyzing {args.input}...")

    # Get the graph statistics and write to a separate file
    num_vertices, num_edges = get_graph_statistics(args.input)
    with open(f"{args.output}/{curr_time}/graph_stats.txt", 'w') as f:
        f.write(f"Name: {args.input}\n")
        f.write(f"Vertices: {num_vertices}\n")
        f.write(f"Edges: {num_edges}\n")
        f.write(f"Unweighted: {args.unweighted}\n")

    for a in algorithms:
        print(f"===============\n{a}\n===============")
        for N in test_sizes:
            output_file = f'{args.output}/{curr_time}/{a}_{N}.txt'

            result = run_cpp_executable(
                N,
                executable=args.executable,
                input_file=args.input,
                output_file=output_file,
                unweighted=args.unweighted,
                algorithm=a
            )

            # Logging output
            print(
                f"Test with N = {N}: {result.seconds} seconds."
            )
            if not result.mismatches:
                print("  PASS")
            else:
                print(f"  FAIL: {result.mismatches} did not output correct distance.")
                all_passed = False

    print(f"Runs saved under {args.output}/{curr_time}.")
    if not all_passed:
        sys.exit(1)


if __name__ == '__main__':
    main()
