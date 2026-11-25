"""
Code identity written into the files this repository produces.

Every file carries the commit that produced it, so the commit must exist
outside the machine that ran it and must describe the code that actually
ran. code_identity() therefore refuses, unless told otherwise, modified or
untracked code and a commit that no remote branch contains.
"""

import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY_URL = "https://github.com/ogourgue/krico-post-production"
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _git(*args):
    try:
        done = subprocess.run(["git", "-C", str(REPOSITORY_ROOT), *args],
                              check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError("git is not available on this node") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"git {' '.join(args)} failed: "
                           f"{exc.stderr.strip()}") from exc
    return done.stdout.strip()


def code_identity(script=None, allow_dirty=False):
    """
    Return the attributes that identify the code: code_repository,
    code_commit (full hash) and code_version (git describe --tags --dirty).

    script is the path of the calling script (__file__). Raises RuntimeError
    when the commit does not describe the code that runs, or does not exist
    on GitHub:
      - tracked files are modified;
      - the calling script, or a file of the krico_recruitment package, is
        not tracked;
      - the commit is on no remote branch (not pushed).
    With allow_dirty=True, for tests, these are accepted and the version
    ends in "-dirty", "-untracked" or "-unpublished".
    """
    commit = _git("rev-parse", "HEAD")
    version = _git("describe", "--tags", "--always", "--dirty")
    problems = []
    modified = _git("status", "--porcelain", "--untracked-files=no")
    if modified:
        problems.append("tracked files are modified:\n" + modified)
    paths = ["recruitment/krico_recruitment"]
    if script is not None:
        paths.append(str(Path(script).resolve().relative_to(REPOSITORY_ROOT)))
    untracked = _git("ls-files", "--others", "--exclude-standard", "--", *paths)
    if untracked:
        problems.append("code files are not tracked:\n" + untracked)
        version += "-untracked"
    if not _git("branch", "--remotes", "--contains", "HEAD"):
        problems.append(f"commit {commit[:7]} is on no remote branch (not pushed)")
        version += "-unpublished"
    if problems and not allow_dirty:
        raise RuntimeError(
            "the commit written in the output would not describe the code "
            "that runs; commit and push first.\n" + "\n".join(problems))
    return {
        "code_repository": REPOSITORY_URL,
        "code_commit": commit,
        "code_version": version,
    }


def history_line(command):
    """One CF history entry: UTC time, then the command that was run."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return f"{now}: {command}"
