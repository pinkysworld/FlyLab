#!/usr/bin/env python3
"""Fail if reachable Git history attributes commits to prohibited AI identities.

FlyLab permits AI-assisted work. This policy only controls *Git attribution* so
repository contributor statistics remain attached to the repository owner.

The check intentionally targets attribution metadata, not ordinary prose. A
commit message may discuss Claude/Anthropic in its body; what is prohibited is
using that identity as author/committer or adding attribution/session trailers.
"""

from __future__ import annotations

import re
import subprocess
import sys
from dataclasses import dataclass


PROHIBITED_NAME = re.compile(r"\bclaude\b", re.IGNORECASE)
PROHIBITED_EMAIL = re.compile(r"@anthropic\.com\b", re.IGNORECASE)
COAUTHOR_TRAILER = re.compile(
    r"^\s*co-authored-by:\s*.*(?:\bclaude\b|@anthropic\.com\b)",
    re.IGNORECASE | re.MULTILINE,
)
SESSION_TRAILER = re.compile(r"^\s*claude-session\s*:", re.IGNORECASE | re.MULTILINE)


@dataclass
class Commit:
    sha: str
    author_name: str
    author_email: str
    committer_name: str
    committer_email: str
    message: str


def git_log() -> list[Commit]:
    # ASCII record/unit separators are safe for Git metadata and make multiline
    # commit messages unambiguous.
    fmt = "%H%x1f%an%x1f%ae%x1f%cn%x1f%ce%x1f%B%x1e"
    out = subprocess.check_output(
        ["git", "log", "--all", f"--format={fmt}"],
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    commits: list[Commit] = []
    for record in out.split("\x1e"):
        record = record.strip("\n")
        if not record.strip():
            continue
        parts = record.split("\x1f", 5)
        if len(parts) != 6:
            raise RuntimeError("Could not parse git log record")
        commits.append(Commit(*parts))
    return commits


def violations(commit: Commit) -> list[str]:
    problems: list[str] = []

    for role, name, email in (
        ("author", commit.author_name, commit.author_email),
        ("committer", commit.committer_name, commit.committer_email),
    ):
        if PROHIBITED_NAME.search(name) or PROHIBITED_EMAIL.search(email):
            problems.append(f"{role}={name} <{email}>")

    if COAUTHOR_TRAILER.search(commit.message):
        problems.append("prohibited Co-Authored-By attribution trailer")
    if SESSION_TRAILER.search(commit.message):
        problems.append("prohibited Claude-Session trailer")

    return problems


def main() -> int:
    bad: list[tuple[Commit, list[str]]] = []

    for commit in git_log():
        problems = violations(commit)
        if problems:
            bad.append((commit, problems))

    if not bad:
        print("Commit attribution policy: OK")
        return 0

    print(
        "Commit attribution policy FAILED.\n"
        "AI-assisted contributions are allowed, but commits must not use "
        "Claude/Anthropic as author, committer, co-author, or session attribution.\n"
        "Rewrite the offending commit metadata while preserving its file tree.\n",
        file=sys.stderr,
    )

    for commit, problems in bad:
        subject = commit.message.splitlines()[0] if commit.message.splitlines() else ""
        print(f"{commit.sha}  {subject}", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
