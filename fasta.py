"""Validated FASTA records for 2-bit sequence archives."""

from dataclasses import dataclass
import re


_INVALID_BASE = re.compile(r"[^ACGTacgt]")
_TOPOLOGY = re.compile(r"(?<![\w-])(circular|linear)(?![\w-])", re.IGNORECASE)
_EXPLICIT_TOPOLOGY = re.compile(r"\btopology\s*=\s*(circular|linear)\b", re.IGNORECASE)


@dataclass(frozen=True)
class FastaRecord:
    name: str
    title: str
    sequence: str
    topology: str = "unknown"


def extract_id(title):
    tokens = title.lstrip('>').split()
    if not tokens:
        return None
    legacy_id = re.search(r'\|([^|]+)\|', tokens[0])
    return legacy_id.group(1) if legacy_id else tokens[0]


def extract_topology(title):
    description = title.split(maxsplit=1)
    if len(description) < 2:
        return "unknown"
    explicit = _EXPLICIT_TOPOLOGY.findall(description[1])
    if explicit:
        values = {value.lower() for value in explicit}
    else:
        text = re.sub(r"\bnon[- ]circular\b", "linear", description[1], flags=re.IGNORECASE)
        values = {value.lower() for value in _TOPOLOGY.findall(text)}
    return values.pop() if len(values) == 1 else "unknown"


def read_fasta(path):
    seen = set()
    name = None
    parts = []
    offset = 0
    with open(path, encoding="utf-8") as file:
        for line_number, line in enumerate(file, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if not parts:
                        raise ValueError(f"{path}:{header_line}: empty FASTA record {name!r}")
                    yield FastaRecord(name, title, ''.join(parts), extract_topology(title))
                title = line[1:].strip()
                name = extract_id(title)
                if name is None:
                    raise ValueError(f"{path}:{line_number}: missing FASTA identifier")
                if '\x00' in title:
                    raise ValueError(f"{path}:{line_number}: NUL byte in FASTA title")
                if name in seen:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA identifier {name!r}")
                seen.add(name)
                header_line = line_number
                parts = []
                offset = 0
            else:
                if name is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
                bad = _INVALID_BASE.search(line)
                if bad:
                    raise ValueError(
                        f"{path}:{line_number}: record {name!r}, base {offset + bad.start() + 1}: "
                        f"unsupported nucleotide {bad.group()!r}; expected A/C/G/T"
                    )
                parts.append(line.upper())
                offset += len(line)
    if name is None:
        raise ValueError(f"{path}: no FASTA records")
    if not parts:
        raise ValueError(f"{path}:{header_line}: empty FASTA record {name!r}")
    yield FastaRecord(name, title, ''.join(parts), extract_topology(title))
