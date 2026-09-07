"""test_readme_links.py - check README's local links and its line anchors

input:   README.md and the files it points at
output:  pytest results
usage:   uv run pytest tests/test_readme_links.py
env:     none; no subprocess calls and no network
notes:   an anchored line must define the name in the link text. a range
         anchor such as #L12-L20 parses as its first line. the end is
         dropped, and README carries no such anchor today
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# [text](target), where target may carry a #L<n> line anchor
LINK_RE = re.compile(r"\[`?([^\]`]+)`?\]\(([^)]+)\)")


def problems(readme):
    """Return one message per broken link, empty when every link is current.
    Paths resolve against the directory holding the README."""
    text = readme.read_text()
    links = []
    for name, url in LINK_RE.findall(text):
        target, _, frag = url.partition("#")
        lineno = int(frag[1:]) if re.fullmatch(r"L\d+", frag) else None
        links.append((name, target, lineno))

    found = []
    # a partial regex failure drops anchors one at a time and would otherwise
    # leave the checks below passing on a shrinking set
    if text.count("#L") != sum(1 for _, _, n in links if n):
        found.append("LINK_RE missed an #L anchor; check it against README")
    for name, target, lineno in links:
        if target.startswith(("http://", "https://")) or not target:
            continue
        path = readme.parent / target
        lines = path.read_text().splitlines() if path.is_file() else None
        if lines is None:
            found.append(f"{target}: not here")
        elif lineno is None:
            continue
        elif lineno > len(lines):
            found.append(f"{target}#L{lineno}: past the end, file has {len(lines)}")
        elif not re.match(rf"\s*(?:(?:def|class)\s+)?{re.escape(name)}\b\s*[(:=]",
                          lines[lineno - 1]):
            found.append(f"{target}#L{lineno}: does not define {name}, "
                         f"{lines[lineno - 1].strip()!r}")
    return found


def test_readme_links_are_current():
    assert not problems(REPO / "README.md")


def test_the_checker_catches_drift(tmp_path):
    (tmp_path / "m.py").write_text("def f():\n    pass\n")
    (tmp_path / "README.md").write_text(
        "[`f`](m.py#L1) [`f`](m.py#L2) [`f`](m.py#L99) [`f`](gone.py) [x](https://a/b)\n"
    )
    assert len(problems(tmp_path / "README.md")) == 3
