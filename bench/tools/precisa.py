"""PRECiSA's bounds. PRECISA names its executable (default ~/.cabal/bin/precisa), patched as the README
says. Branch and bound runs deeper than the defaults (-d 7 -p 14), which leave relative bounds loose.
It takes binary64 only. Some runs grow by gigabytes a second, so a run stops past MEMORY bytes and
fewer run at once."""

import ctypes
import decimal
import json
import os
import signal
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

from milo.fpcore import fpcore, parse

from .common import cached

PRECISA = os.path.expanduser(os.environ.get("PRECISA", "~/.cabal/bin/precisa"))
OPTIONS = ["--fpcore", "--relative-error", "--json", "-d", "9", "-p", "16"]
MEMORY, JOBS, CHUNK = 2 << 30, 4, 32
LIBC = ctypes.CDLL("/usr/lib/libSystem.B.dylib")


class RusageInfo(ctypes.Structure):
    """macOS's rusage_info_v0, filled by proc_pid_rusage."""
    _fields_ = [("uuid", ctypes.c_uint8 * 16)] + [(f, ctypes.c_uint64) for f in (
        "user_time", "system_time", "pkg_idle_wkups", "interrupt_wkups", "pageins", "wired_size",
        "resident_size", "phys_footprint", "proc_start_abstime", "proc_exit_abstime")]


def footprint(pid):
    """Bytes of memory a process uses, counting compressed and swapped pages, which its resident size
    leaves out."""
    info = RusageInfo()
    return info.phys_footprint if LIBC.proc_pid_rusage(pid, 0, ctypes.byref(info)) == 0 else 0


def run(cmd, timeout, cwd):
    """(output, seconds, status) of cmd, stopped at the timeout or past MEMORY bytes."""
    start, status = time.time(), "ok"
    with tempfile.TemporaryFile("w+") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, text=True, cwd=cwd, start_new_session=True)
        while proc.poll() is None:
            status = "timeout" if time.time() - start > timeout else "memory" if footprint(proc.pid) > MEMORY else "ok"
            if status != "ok":
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except OSError:  # it exited meanwhile
                    pass
                proc.wait()
                break
            time.sleep(0.05)
        log.seek(0)
        return log.read(), round(time.time() - start, 2), status


def bound_of(x):
    """A bound from PRECiSA's report, which writes a missing one as "infinity" or null."""
    return float(x) if isinstance(x, (int, float)) else None


def renamed(core):
    """The core with its inputs named v0, v1, ..., since PRECiSA's lexer rejects names such as x.re."""
    names = {x: f"v{i}" for i, x in enumerate(core.args)}

    def rename(e):
        match e:
            case ("Var", x):
                return ("Var", names[x])
            case ("Num", _):
                return e
            case (op, *args):
                return (op, *map(rename, args))

    return core._replace(args=list(names.values()), box={names[x]: b for x, b in core.box.items()},
                         body=rename(core.body))


def one(text, workdir, timeout):
    core = renamed(parse(text))
    if core.precision != "binary64":
        return {"status": "unsupported", "abs": None, "rel": None}
    d = tempfile.mkdtemp(dir=workdir)  # PRECiSA writes its certificates next to its input
    with open(os.path.join(d, "p.fpcore"), "w") as f:
        f.write(fpcore(core).replace("(FPCore (", "(FPCore p (", 1))
    with open(os.path.join(d, "p.input"), "w") as f:  # ranges in decimal, since it takes no exponents
        ranges = (f"{x} in [{decimal.Decimal(lo):f}, {decimal.Decimal(hi):f}]" for x, (lo, hi) in core.box.items())
        f.write(f"p({','.join(core.args)}): {', '.join(ranges)}")
    out, seconds, status = run([PRECISA, "p.fpcore", "p.input", *OPTIONS], timeout, d)
    report = {"seconds": seconds, "status": status, "abs": None, "rel": None}
    if status == "ok" and '{"results"' in out:
        result = json.loads(out[out.index('{"results"'):])["results"][0]
        if result.get("unstableError") is None:  # a straight-line program has only the stable path
            report["abs"] = bound_of(result.get("stableError"))
            report["rel"] = bound_of(result.get("relativeStableError"))
    if status == "ok" and report["abs"] is None:
        report["status"] = "fail"
    return report


def bound(texts, workdir, timeout, jobs):
    """Each text's report, cached in chunks so that an interrupted run keeps what it finished."""
    def work(todo):
        with ThreadPoolExecutor(min(jobs, JOBS)) as pool:
            return dict(zip(todo, pool.map(lambda t: one(t, workdir, timeout), todo)))

    out = {}
    for i in range(0, len(texts), CHUNK):
        out.update(cached(f"precisa exact-inputs exact-constants {OPTIONS} {timeout}", texts[i:i + CHUNK], work))
    return out
