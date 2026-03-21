from __future__ import annotations

import ssl
import tarfile
from pathlib import Path
from urllib.request import Request, urlopen


PBMC3K_URL = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
PBMC3K_ARCHIVE = "pbmc3k_filtered_gene_bc_matrices.tar.gz"


def download_with_ssl_fallback(url: str, destination: Path) -> None:
    request = Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/122.0.0.0 Safari/537.36"
            )
        },
    )
    try:
        with urlopen(request) as response, destination.open("wb") as handle:
            handle.write(response.read())
    except ssl.SSLCertVerificationError:
        insecure_context = ssl._create_unverified_context()
        with urlopen(request, context=insecure_context) as response, destination.open("wb") as handle:
            handle.write(response.read())
    except Exception as exc:
        message = str(exc).lower()
        if "certificate verify failed" not in message:
            raise
        insecure_context = ssl._create_unverified_context()
        with urlopen(request, context=insecure_context) as response, destination.open("wb") as handle:
            handle.write(response.read())


def download_dataset(data_dir: Path) -> Path:
    data_dir.mkdir(parents=True, exist_ok=True)
    archive_path = data_dir / PBMC3K_ARCHIVE
    matrix_dir = find_10x_matrix_dir(data_dir)

    if matrix_dir is not None:
        return matrix_dir

    if not archive_path.exists():
        print("Downloading PBMC 3k single-cell dataset from 10x Genomics...")
        download_with_ssl_fallback(PBMC3K_URL, archive_path)

    print("Extracting dataset archive...")
    with tarfile.open(archive_path, "r:gz") as handle:
        handle.extractall(path=data_dir)

    matrix_dir = find_10x_matrix_dir(data_dir)
    if matrix_dir is None:
        raise FileNotFoundError(
            f"Could not locate an extracted 10x matrix directory under {data_dir}"
        )

    return matrix_dir


def find_10x_matrix_dir(root: Path) -> Path | None:
    candidates = []
    for matrix_file in root.rglob("matrix.mtx"):
        parent = matrix_file.parent
        has_barcodes = (parent / "barcodes.tsv").exists() or (parent / "barcodes.tsv.gz").exists()
        has_features = (
            (parent / "genes.tsv").exists()
            or (parent / "genes.tsv.gz").exists()
            or (parent / "features.tsv").exists()
            or (parent / "features.tsv.gz").exists()
        )
        if has_barcodes and has_features:
            candidates.append(parent)

    if not candidates:
        return None

    candidates = sorted(candidates, key=lambda path: len(path.parts))
    return candidates[0]
