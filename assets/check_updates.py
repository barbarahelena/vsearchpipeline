#!/usr/bin/env python3
"""
check_updates.py — check Bioconda/conda-forge packages for new versions
and send an email notification if any updates are found.

Environment variables required (only when updates are found):
  GMAIL_USER         — sender Gmail address
  GMAIL_APP_PASSWORD — Gmail App Password (not your account password)
  NOTIFY_EMAIL       — recipient address

Usage:
  python3 check_updates.py [--tools tools.txt] [--versions versions_seen.json]
"""

import argparse
import json
import os
import smtplib
import sys
from email.mime.text import MIMEText
from pathlib import Path

import requests

CHANNELS = ["bioconda", "conda-forge"]
ANACONDA_API = "https://api.anaconda.org/package/{channel}/{package}"


def get_latest_version(package: str) -> tuple[str | None, str | None]:
    """Return (channel, version) for the latest release, trying channels in order."""
    for channel in CHANNELS:
        url = ANACONDA_API.format(channel=channel, package=package)
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 404:
                continue
            r.raise_for_status()
            data = r.json()
            version = data.get("latest_version") or data.get("version")
            if version:
                return channel, version
        except requests.RequestException as exc:
            print(f"  WARNING: HTTP error for {package} on {channel}: {exc}", file=sys.stderr)
            continue
    return None, None


def send_email(updates: list[dict], gmail_user: str, app_password: str, notify_email: str) -> None:
    lines = ["The following pipeline tool updates are available:\n"]
    for u in updates:
        old = u["old"] or "(not tracked)"
        lines.append(f"  {u['package']:35s}  {old:15s} -> {u['new']}  [{u['channel']}]")
    lines.append("\nUpdate the version pins in conf/modules.config and the Nextflow process definitions accordingly.")
    body = "\n".join(lines)

    msg = MIMEText(body)
    msg["Subject"] = f"[vsearchpipeline] {len(updates)} tool update(s) available"
    msg["From"] = gmail_user
    msg["To"] = notify_email

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as smtp:
        smtp.login(gmail_user, app_password)
        smtp.sendmail(gmail_user, notify_email, msg.as_string())
    print(f"Email sent to {notify_email}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tools",    default="assets/tools.txt",          help="Path to tools list")
    parser.add_argument("--versions", default="assets/versions_seen.json", help="Path to versions_seen.json")
    args = parser.parse_args()

    tools_path    = Path(args.tools)
    versions_path = Path(args.versions)

    # Load tool list
    packages = [
        line.strip()
        for line in tools_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    print(f"Checking {len(packages)} packages...")

    # Load known versions
    if versions_path.exists():
        known = json.loads(versions_path.read_text())
    else:
        known = {}

    updates = []
    new_known = dict(known)

    for package in packages:
        channel, latest = get_latest_version(package)
        if latest is None:
            print(f"  WARNING: {package} not found on any channel — skipping")
            continue

        old = known.get(package)
        print(f"  {package}: known={old or '(new)'}  latest={latest}  [{channel}]")

        if latest != old:
            updates.append({"package": package, "old": old, "new": latest, "channel": channel})

        new_known[package] = latest

    # Write updated versions regardless of whether there were updates
    versions_path.write_text(json.dumps(new_known, indent=2) + "\n")
    print(f"\nWrote updated {versions_path}")

    if not updates:
        print("No updates found.")
        return

    print(f"\n{len(updates)} update(s) found:")
    for u in updates:
        print(f"  {u['package']}: {u['old'] or '(new)'} -> {u['new']}  [{u['channel']}]")

    # Send email notification
    gmail_user   = os.environ.get("GMAIL_USER")
    app_password = os.environ.get("GMAIL_APP_PASSWORD")
    notify_email = os.environ.get("NOTIFY_EMAIL")

    if not all([gmail_user, app_password, notify_email]):
        print("\nWARNING: Email credentials not set — skipping notification.")
        print("Set GMAIL_USER, GMAIL_APP_PASSWORD, and NOTIFY_EMAIL to enable email.")
        return

    send_email(updates, gmail_user, app_password, notify_email)


if __name__ == "__main__":
    main()
