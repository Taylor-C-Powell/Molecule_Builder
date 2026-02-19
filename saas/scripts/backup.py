"""SQLite database backup utility for MolBuilder SaaS.

Uses the SQLite Online Backup API (sqlite3.backup) for safe, consistent
backups even while the application is running with WAL mode.

Usage:
    # Backup to local directory
    python -m scripts.backup

    # Backup to specific directory
    python -m scripts.backup --output /path/to/backups

    # Backup with retention (delete backups older than N days)
    python -m scripts.backup --retain-days 30
"""

import argparse
import logging
import os
import sqlite3
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("molbuilder.backup")

# Default database paths (match app/config.py defaults)
DEFAULT_DBS = [
    "molbuilder_users.db",
    "molbuilder_audit.db",
    "molbuilder_molecules.db",
]


def backup_database(src_path: str, dst_path: str) -> bool:
    """Backup a SQLite database using the Online Backup API.

    This is safe to call while the database is being written to.
    The backup API handles WAL mode correctly.
    """
    if not os.path.exists(src_path):
        logger.warning("Source database not found: %s", src_path)
        return False

    try:
        src = sqlite3.connect(src_path)
        dst = sqlite3.connect(dst_path)
        with dst:
            src.backup(dst)
        dst.close()
        src.close()

        src_size = os.path.getsize(src_path)
        dst_size = os.path.getsize(dst_path)
        logger.info(
            "Backed up %s -> %s (%d KB)",
            src_path, dst_path, dst_size // 1024,
        )

        # Verify the backup is readable
        verify = sqlite3.connect(dst_path)
        verify.execute("PRAGMA integrity_check")
        verify.close()

        return True
    except Exception:
        logger.exception("Failed to backup %s", src_path)
        return False


def cleanup_old_backups(backup_dir: str, retain_days: int) -> int:
    """Delete backup files older than retain_days."""
    cutoff = datetime.now(timezone.utc) - timedelta(days=retain_days)
    removed = 0
    for f in Path(backup_dir).glob("*.db"):
        if f.stat().st_mtime < cutoff.timestamp():
            f.unlink()
            removed += 1
            logger.info("Removed old backup: %s", f)
    return removed


def run_backup(output_dir: str, db_dir: str, retain_days: int | None = None) -> bool:
    """Run backup for all MolBuilder databases."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    os.makedirs(output_dir, exist_ok=True)

    success = True
    for db_name in DEFAULT_DBS:
        src = os.path.join(db_dir, db_name)
        dst = os.path.join(output_dir, f"{timestamp}_{db_name}")
        if not backup_database(src, dst):
            success = False

    if retain_days is not None:
        removed = cleanup_old_backups(output_dir, retain_days)
        if removed:
            logger.info("Cleaned up %d old backup(s)", removed)

    return success


def main():
    parser = argparse.ArgumentParser(description="Backup MolBuilder SQLite databases")
    parser.add_argument(
        "--output", "-o",
        default=os.environ.get("BACKUP_DIR", "./backups"),
        help="Backup output directory (default: ./backups or $BACKUP_DIR)",
    )
    parser.add_argument(
        "--db-dir",
        default=os.environ.get("DB_DIR", "."),
        help="Directory containing source databases (default: . or $DB_DIR)",
    )
    parser.add_argument(
        "--retain-days",
        type=int,
        default=None,
        help="Delete backups older than N days",
    )
    args = parser.parse_args()

    logger.info("Starting MolBuilder database backup")
    ok = run_backup(args.output, args.db_dir, args.retain_days)
    if ok:
        logger.info("Backup completed successfully")
    else:
        logger.error("Backup completed with errors")
        sys.exit(1)


if __name__ == "__main__":
    main()
