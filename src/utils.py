# =============================================================
# Shared utilities — config loader, logger, path helpers
# Import this in every script: from src.utils import load_config, get_logger
# =============================================================

import yaml
from pathlib import Path
from loguru import logger
import sys


# Project Root
# Always resolves correctly regardless of where script is run from
PROJECT_ROOT = Path(__file__).parent.parent


# Config Loader
def load_config() -> dict:
    """
    Load and return the central settings.yaml configuration.
    Returns a nested dict matching the YAML structure.
    """
    config_path = PROJECT_ROOT / "configs" / "settings.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found at: {config_path}")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# Logger Setup
def get_logger(name: str):
    """
    Returns a configured loguru logger.
    Logs to both console and file (path set in settings.yaml).

    Usage:
        from src.utils import get_logger
        log = get_logger(__name__)
        log.info("Starting ingestion...")
    """
    cfg = load_config()
    log_cfg = cfg.get("logging", {})

    # Remove default logger
    logger.remove()

    # Console logger — coloured, readable
    logger.add(
        sys.stdout,
        level=log_cfg.get("level", "INFO"),
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
               "<level>{level: <8}</level> | "
               "<cyan>{name}</cyan> | "
               "<level>{message}</level>",
        colorize=True
    )

    # File logger — plain text, rotated
    if log_cfg.get("log_to_file", False):
        log_file = PROJECT_ROOT / log_cfg.get("log_file", "outputs/logs/pipeline.log")
        log_file.parent.mkdir(parents=True, exist_ok=True)
        logger.add(
            log_file,
            level=log_cfg.get("level", "INFO"),
            rotation=log_cfg.get("rotation", "10 MB"),
            retention=log_cfg.get("retention", "7 days"),
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name} | {message}"
        )

    return logger.bind(name=name)


# Path Helpers
def get_path(key_path: str) -> Path:
    """
    Resolve a dot-notation path key from settings.yaml to an absolute Path.

    Examples:
        get_path("paths.raw.nomenclature.depmap_sample")
        get_path("paths.processed.master_index")
    """
    cfg = load_config()
    keys = key_path.split(".")
    value = cfg
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise KeyError(f"Key '{key}' not found in config path: '{key_path}'")
        value = value[key]
    return PROJECT_ROOT / value


def ensure_processed_dir():
    """
    Create data/processed/ and outputs/ directories if they don't exist.
    Call this at the start of any script that writes output files.
    """
    cfg = load_config()
    dirs = [
        PROJECT_ROOT / cfg["paths"]["processed"]["base"],
        PROJECT_ROOT / cfg["paths"]["outputs"]["base"],
        PROJECT_ROOT / cfg["paths"]["outputs"]["reports"],
        PROJECT_ROOT / cfg["paths"]["outputs"]["exports"],
        PROJECT_ROOT / cfg["paths"]["outputs"]["logs"],
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


# DataFrame Helpers
def save_parquet(df, key_path: str, log=None):
    """
    Save a DataFrame to parquet using path from settings.yaml.

    Usage:
        save_parquet(df, "paths.processed.master_index")
    """
    out_path = get_path(key_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=True)
    if log:
        log.info(f"Saved {len(df):,} rows → {out_path.relative_to(PROJECT_ROOT)}")


def load_parquet(key_path: str):
    """
    Load a parquet file using path from settings.yaml.

    Usage:
        df = load_parquet("paths.processed.master_index")
    """
    import pandas as pd
    in_path = get_path(key_path)
    if not in_path.exists():
        raise FileNotFoundError(
            f"Processed file not found: {in_path}\n"
            f"Have you run the pipeline up to this stage?"
        )
    return pd.read_parquet(in_path)


# Quick Sanity Check 
if __name__ == "__main__":
    # Run this directly to verify everything is wired up correctly:
    # python -m src.utils
    log = get_logger(__name__)
    cfg = load_config()

    log.info(f"Project root: {PROJECT_ROOT}")
    log.info(f"Config loaded: {cfg['project']['name']} v{cfg['project']['version']}")

    ensure_processed_dir()
    log.info("Processed and output directories verified.")

    # Test path resolution
    test_path = get_path("paths.raw.nomenclature.depmap_sample")
    log.info(f"Sample path resolves to: {test_path}")
    log.info(f"File exists: {test_path.exists()}")

    log.success("utils.py — all checks passed ✓")
