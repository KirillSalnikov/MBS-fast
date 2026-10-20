"""Read sampling metadata from legacy and compact fullauto checkpoints."""
from pathlib import Path


def fnv(data):
    value = 14695981039346656037
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return format(value, "x")


def sampling_text(folder):
    folder = Path(folder)
    legacy = folder / "sampling.index"
    if legacy.exists():
        data = legacy.read_bytes()
        checksum = folder / "sampling.index.fnv"
        if checksum.exists() and fnv(data) != checksum.read_text().strip():
            raise ValueError("sampling index checksum mismatch: " + str(folder))
        return data.decode("utf-8")
    marker = (folder / "complete.txt").read_bytes()
    try:
        _, header, data = marker.split(b"\n", 2)
        magic, checksum, length = header.split()
        size = int(length)
    except (ValueError, TypeError) as error:
        raise ValueError("invalid compact sampling header: " + str(folder)) from error
    if magic != b"MBS_SAMPLING_V1" or size < 0 or size > len(data) or checksum != fnv(data[:size]).encode("ascii"):
        raise ValueError("compact sampling checksum/size mismatch: " + str(folder))
    tail = data[size:]
    if tail:
        try:
            header, diagnostics = tail.split(b"\n", 1)
            magic, checksum = header.split()
        except ValueError as error:
            raise ValueError("invalid compact diagnostics header: " + str(folder)) from error
        if magic != b"MBS_DIAGNOSTICS_V1" or checksum != fnv(diagnostics).encode("ascii"):
            raise ValueError("compact diagnostics checksum mismatch: " + str(folder))
    return data[:size].decode("utf-8")
