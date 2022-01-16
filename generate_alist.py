"""
generate_alist.py

Creates .alist files for the LDPC codes outlined in the IRIG106 Telemetry Standards
The generated files are quite large, and do not exist within the gnuradio/gnuradio repo
on github.  Since I wanted to experiment with these particular LDPC codes, I've made this
script to generate them.

I couldn't be bothered hitting ctrl+C and ctrl+V thousands of times so this script 
downloads the relevant PDF file, grabs the information required for populating generator 
matricies, parses said data, creates corresponding generator matricies, and finally writes 
them to file in the alist format (http://www.inference.org.uk/mackay/codes/alist.html).

Author: C S McDiarmid
"""

import logging
from typing import List, Dict, Tuple
from pathlib import Path
import re

import numpy as np
import pandas as pd
import tabula
import requests


IRIG106_2020_CH2_URL = "http://www.irig106.org/docs/106-20/chapter2.pdf"
RELEVANT_CH2_PAGES = list(range(107, 138))
TABLE_TITLE_IRIG106_PATTERN = re.compile(r"""
    Table D-\d. ?
    First Rows of Circulants in Generator Matrix
    , r=\d+/\d+, k=\d+
    """,
    re.VERBOSE
)
PARAM_IRIG106_PATTERN = re.compile(r'r=(?P<r>\d+/\d+), k=(?P<k>\d+)')
REPO_DIR = Path(__file__).parent

_logger = logging.getLogger(__name__)


def ishexadecimal(string: str) -> bool:
    return bool(re.fullmatch(r'([0-9a-fA-F]+)', string))


def download_pdf(
        url: str,
        destination: Path
) -> Path:
    """
    Downloads a requested pdf, names the file based on the tail of the url, 
    and saves to a provided destination directory

    :param url: URL of the desired PDF file
    :return: Path object representing the location of the downloaded PDF
    """
    # Nab the filename from the url
    *_, filename = url.split('/')
    if not filename.endswith('.pdf'):
        filename += '.pdf'
    filepath = destination / filename

    if filepath.exists():
        _logger.info(f"{filepath} already exists, skipping download.")
        return filepath

    # Request the file
    _logger.info(f"Requesting {url}")
    response = requests.get(url)
    response.raise_for_status()
    
    # Ensure we can write this file to the specified directory
    destination.mkdir(parents=True, exist_ok=True)
    bytes_written = filepath.write_bytes(response.content)
    _logger.info(f"Wrote {bytes_written} bytes to {filepath}")
    return filepath


def construct_matrix(
        k: int,
        n: int,
        circulant_size: int,
        circulants: List[int]
) -> np.array:
    """
    Constructs an LDPC generator matrix based off size info and a list of circulants

    :param k: Generator matrix rows
    :param n: Generator matrix columns
    :param circulant_size: Circulant matrix size (number of bits in each circulant)
    :param circulants: List of circulants in integer format
    :return: Generator matrix
    """
    # Allocate the generator matrix and calculate
    _logger.info(
        f"Constructing generator matrix of size k={k} n={n} "
        f"containing {len(circulants)} circulant sub-matricies "
        f"of size {circulant_size}"
    )
    matrix = np.zeros((k, n))
    matrix[:k, :k] = np.identity(k)
    columns_of_circulants = (n - k) // circulant_size
    rows_of_circulants = k // circulant_size

    for i, circ_row in enumerate(circulants):
        # Chunk row and column
        sub_row, sub_col = divmod(i, columns_of_circulants)
        _logger.debug(f"Creating circulant matrix ({sub_row}, {sub_col})")

        # Array used to generate each row of this chunk's matrix
        row_array = np.array([int(x) for x in f'{circ_row:0{circulant_size}b}'])
        for row in range(circulant_size):
            matrix[
                row + sub_row * circulant_size,
                k + sub_col * circulant_size:k + (sub_col+1) * circulant_size
            ] = np.roll(row_array, row)

    return matrix


def hex_table_to_circulant_rows(table: List[str]) -> Tuple[int, List[int]]:
    """
    Converts circulant rows in string hex format to integer format

    :param table: List of circulants in string hex format
    :return: Circulant submatrix size, circulant row 1 integer format
    """
    # Each hex digit = 4 bits
    if not all(len(circ) == len(table[0]) and ishexadecimal(circ) for circ in table):
        raise ValueError("All items must be of equal length and hexadecimal format.")
    
    _logger.info("Validated table of hexadecimal circulant rows.")
    return len(table[0]) * 4, [int(circ, 16) for circ in table]


def extract_generator_matricies_irig106() -> Dict[str, np.array]:
    """
    Downloads IRIG106 Chapter 2, extracts information pertaining to LDPC
    generator matricies, constructs generator matricies and returns them

    :return: LDPC generator matricies indexed by name
    """
    # Download the PDF file (if it doesn't already exist)
    pdf_file = download_pdf(
        url=IRIG106_2020_CH2_URL,
        destination=REPO_DIR / "documents"
    )
    
    # Extract relevant tabulated data
    dfs: List[pd.DataFrame] = tabula.read_pdf(
        pdf_file,
        pages=RELEVANT_CH2_PAGES,
        lattice=True  # Required for one-row tables continued from prev page
    )

    # Join relevant page separated dfs pertaining to the same table
    generator_dfs: Dict[str, List[str]] = {}
    key = None

    for i, df in enumerate(dfs):
        if len(df.columns) < 2 and df.empty:  # Fix for random empty tables
            continue

        if df.empty:  # Fix for single row tables (from previous page)
            key_column = df.columns[1]
            unfiltered_data = []
        else:
            key_column = df.columns[-2]
            unfiltered_data = [x.replace('\r', '') for x in df[key_column].values]
        
        if key_column.startswith('Table'):
            key = key_column

            if "Row 1" in unfiltered_data:
                generator_start = unfiltered_data.index("Row 1")
            else:
                generator_start = -1
            
            key += ''.join(f' {item}' for item in unfiltered_data[:generator_start])
            generator_dfs[key] = []
            del unfiltered_data[:generator_start]
        else:
            unfiltered_data.insert(0, key_column.replace('\r', ''))

        if TABLE_TITLE_IRIG106_PATTERN.search(key):
            generator_dfs[key].extend([x for x in unfiltered_data if ishexadecimal(x)])
        else:
            del generator_dfs[key]

    # Create generator matricies from tables of circulant rows
    generators: Dict[str, np.array] = {}
    for title, table in generator_dfs.items():
        match = PARAM_IRIG106_PATTERN.search(title)
        params = match.groupdict()
        k = int(params.get('k'))
        r = [int(x) for x in params.get('r').split('/')]
        n = k * r[1] // r[0]
        name = f"LDPC r={r} n={n} k={k}"
        generators[name] = construct_matrix(k, n, *hex_table_to_circulant_rows(table))

    return generators


def write_matrix_to_alist(
        name: str,
        matrix: np.array
) -> None:
    """
    Writes an alist file for a given LDPC generator matrix

    :param name: Name/title of the alist file
    :param matrix: LDPC generator matrix
    """
    # Iteration information
    k, n = matrix.shape
    col_sums = [sum(matrix[:, col]) for col in range(n)]
    row_sums = [sum(matrix[row,:]) for row in range(k)]

    # Prepare data in alist format: http://www.inference.org.uk/mackay/codes/alist.html
    ordered_iterables = [
        matrix.shape[::-1],
        (max(col_sums), max(row_sums)),
        col_sums,
        row_sums,
        *[np.where(matrix[:, col])[0] for col in range(n)],
        *[np.where(matrix[row, :])[0] for row in range(k)],
    ]

    # Dump lines to file
    name = re.sub(r"[\[\]\(\)/=, ]+", "_", name).lower()
    filepath = REPO_DIR / "alist" / f"{name}.alist"
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        f.writelines([' '.join(f'{x:.0f}' for x in it) +  ' \n' for it in ordered_iterables])
        _logger.info(f"Wrote {filepath.stat().st_size} bytes to {filepath}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    generator_matricies = extract_generator_matricies_irig106()
    for name, matrix in generator_matricies.items():
        write_matrix_to_alist(name, matrix)
