#!/usr/bin/env python3
"""
Automated pre-release checks for GSL.

This runs the checks listed in HACKING ("Pre-release checks" and "Notes on
making a release"), and some more, on a clean export of a git revision
(default HEAD), building everything out of tree in a work directory: the
working tree is never configured, modified or built.

Checks (see --list):

  static      exported symbols prefixed with gsl_/GSL_ (cblas_ for CBLAS);
              config.h use in sources and installed headers; all headers
              compiled together (C, C++) and each header on its own (C89
              -pedantic, C99 with HAVE_INLINE, C11 after <complex.h>, C++);
              version numbers consistent (configure.ac, NEWS, doc/conf.py);
              autoreconf -W all
  build:*     configure, make, make check: strict C89 (which uses all the
              math fallbacks), strict C99, C2x, the fallbacks with a C99
              compiler, -O3 with FMA contraction, 32-bit x86, without
              HAVE_INLINE, with GSL_DISABLE_DEPRECATED, ASan+UBSan, and
              clang/icc if present; implicit declarations are failures
  dynamic     GSL_TEST_VERBOSE=1 make check; example programs
              (scripts/test_examples.py); valgrind memcheck of every test
              program (static build, so that the test programs are not
              libtool wrapper scripts)
  dist        make distcheck (default and with install-sh); make dist,
              unpack, build, and make check under a virtual memory limit
  docs        Sphinx html, latexpdf, info
  abi         compare with the previous release (abidiff and nm): removed,
              added and changed interfaces, and whether the libtool
              version numbers in configure.ac are consistent with them
  qemu:*      cross-compile and run every test program under QEMU user
              mode (aarch64, ppc64le, ppc64, sparc64, ...), for every
              target whose cross compiler and qemu are installed

A full run takes many hours (valgrind and the QEMU runs dominate). Status
lines are printed as checks start and finish; every command's output goes to
a per-check log; a summary (summary.txt, summary.json) is written at the end
and on interruption. --resume skips checks which passed in a previous run in
the same work directory.

Examples:

  scripts/prerelease.py --list
  scripts/prerelease.py --workdir /scratch/gsl-rc --jobs 4 --nice 19
  scripts/prerelease.py --quick
  scripts/prerelease.py --only 'build:c89*,headers' --rev release-2-9
  scripts/prerelease.py --resume --skip valgrind
"""

import argparse
import fnmatch
import glob
import json
import os
import platform
import re
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import time

# ---------------------------------------------------------------------------
# results


class Result:
    """outcome of a check: PASS, WARN (worth a look), FAIL, SKIP"""

    def __init__(self, status, summary="", details=None):
        self.status = status
        self.summary = summary
        self.details = details or []    # lines worth showing in the summary

    def to_dict(self):
        return {"status": self.status, "summary": self.summary, "details": self.details}

    @staticmethod
    def from_dict(d):
        return Result(d["status"], d.get("summary", ""), d.get("details", []))


def worst(*statuses):
    order = ["PASS", "SKIP", "WARN", "FAIL"]
    return max(statuses, key=order.index)


class CheckError(Exception):
    """a command a check depends on failed"""

    def __init__(self, msg, details=None):
        Exception.__init__(self, msg)
        self.details = details or []


# ---------------------------------------------------------------------------
# running commands


class Context:
    def __init__(self, args):
        self.args = args
        self.workdir = os.path.abspath(args.workdir)
        self.jobs = args.jobs
        self.repo = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        os.makedirs(self.workdir, exist_ok=True)
        self.master_log = open(os.path.join(self.workdir, "prerelease.log"), "a")
        self.src = None             # exported source tree
        self.commit = None
        self.builds = {}            # name -> build dir of successful builds
        self.current_log = None

    # status output: terminal and master log
    def status(self, msg):
        line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
        print(line, flush=True)
        self.master_log.write(line + "\n")
        self.master_log.flush()

    def run(self, cmd, cwd, log=None, env=None, timeout=None, check=True, shell=False, quiet=False):
        """run a command, appending its output to log; return the exit code
        (raises CheckError if check and the command fails)"""
        log = log or self.current_log
        timeout = timeout or self.args.timeout * 3600
        e = dict(os.environ)
        e["LC_ALL"] = "C"
        if env:
            e.update(env)
        text = cmd if isinstance(cmd, str) else " ".join(shlex.quote(c) for c in cmd)
        with open(log, "a") as f:
            f.write("\n$ cd %s && %s\n" % (cwd, text))
            if env:
                f.write("  env: %s\n" % " ".join("%s=%s" % kv for kv in env.items()))
            f.flush()
            t0 = time.time()
            try:
                p = subprocess.run(cmd, cwd=cwd, env=e, stdout=f, stderr=subprocess.STDOUT,
                                   timeout=timeout, shell=shell, preexec_fn=self._preexec)
                rc = p.returncode
            except subprocess.TimeoutExpired:
                f.write("*** timeout after %d s\n" % timeout)
                rc = -9999
            f.write("[exit %d, %.0f s]\n" % (rc, time.time() - t0))
        if check and rc != 0:
            raise CheckError("command failed (exit %d): %s" % (rc, text),
                             ["see %s" % log] + tail(log, 15))
        return rc

    def _preexec(self):
        # own process group, so that a timeout or ^C reaches the whole tree
        os.setpgrp()

    def make(self, cwd, targets=(), env=None, jobs=None, check=True, keep_going=False, extra=()):
        cmd = ["make", "-j%d" % (jobs or self.jobs)]
        if keep_going:
            cmd.append("-k")
        cmd += list(extra) + list(targets)
        return self.run(cmd, cwd, env=env, check=check)


def tail(path, n):
    try:
        with open(path, errors="replace") as f:
            return [l.rstrip("\n") for l in f.readlines()[-n:]]
    except OSError:
        return []


def which(prog):
    return shutil.which(prog)


def is_elf(path):
    try:
        with open(path, "rb") as f:
            return f.read(4) == b"\x7fELF"
    except OSError:
        return False


def human(seconds):
    seconds = int(seconds)
    if seconds < 60:
        return "%ds" % seconds
    if seconds < 3600:
        return "%dm%02ds" % (seconds // 60, seconds % 60)
    return "%dh%02dm" % (seconds // 3600, (seconds % 3600) // 60)


# ---------------------------------------------------------------------------
# source tree


def prepare_source(ctx):
    """export the revision to be checked and generate the build system"""
    if ctx.src:
        return ctx.src
    rev = ctx.args.rev
    commit = subprocess.run(["git", "rev-parse", "--short", rev], cwd=ctx.repo,
                            capture_output=True, text=True, check=True).stdout.strip()
    src = os.path.join(ctx.workdir, "src-%s" % commit)
    ok = os.path.join(src, ".prerelease-ok")
    log = os.path.join(ctx.workdir, "source.log")
    if not os.path.exists(ok):
        ctx.status("exporting %s (%s) and running autogen.sh" % (rev, commit))
        if os.path.exists(src):
            shutil.rmtree(src)
        os.makedirs(src)
        ctx.run("git archive %s | tar -x -C %s" % (shlex.quote(rev), shlex.quote(src)),
                ctx.repo, log=log, shell=True)
        ctx.run(["sh", "autogen.sh"], src, log=log)
        open(ok, "w").close()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ctx.repo,
                           capture_output=True, text=True).stdout.strip()
    if dirty and rev == "HEAD":
        ctx.status("note: the working tree has uncommitted changes, which are NOT tested")
    ctx.src, ctx.commit = src, commit
    return src


# ---------------------------------------------------------------------------
# build configurations

IS_X86 = platform.machine() in ("x86_64", "i686", "i386", "AMD64")

# the C99/POSIX math functions which have fallbacks in the library
MATH_FUNCS = ["erf", "erfc", "exp2", "expm1", "fma", "lgamma", "lgamma_r", "log1p", "tgamma"]

# name: (description, configure arguments, config.h edits, run make check?)
#
# Most variations concern how the compiler sees the source (language
# standard, inline functions, deprecated declarations, fallbacks,
# sanitizers): they are independent of the CPU and are run on the build
# machine only. What differs between architectures (endianness, long double
# format, FMA contraction, char signedness, alignment) is covered by one
# default build per architecture under QEMU (the qemu:* checks).
BUILDS = {
    "default":        ("default flags", [], None, True),
    "c89":            ("strict C89: CFLAGS=\"-g -O2 -std=c89 -pedantic -Wall\" (all math fallbacks)",
                       ["CFLAGS=-g -O2 -std=c89 -pedantic -Wall"], None, True),
    "c99":            ("strict C99: CC=c99", ["CC=c99"], None, True),
    "c2x":            ("-std=c2x (new keywords, no implicit declarations)", ["CFLAGS=-g -O2 -std=c2x"], None, True),
    "fallbacks":      ("library fallbacks for the C99 math functions (ac_cv_func_*=no)",
                       ["ac_cv_func_%s=no" % f for f in MATH_FUNCS], None, True),
    "O3-fma":         ("-O3 with FMA contraction: -march=x86-64-v3 (x86)", ["CFLAGS=-g -O3 -march=x86-64-v3"], None, True),
    # the i686 cross compiler rather than gcc -m32: on Debian and Ubuntu,
    # gcc-multilib conflicts with the cross compilers of the qemu:* checks.
    # Linked statically, the tests run natively on an x86-64 kernel
    "x86-32":         ("32-bit x86: CC=\"i686-linux-gnu-gcc -static\", run natively "
                       "(configure selects SSE2 arithmetic instead of x87)",
                       ["--host=i686-linux-gnu", "--disable-shared", "CC=i686-linux-gnu-gcc -static"], None, True),
    "no-inline":      ("HAVE_INLINE undefined in config.h", [], {"HAVE_INLINE": None, "HAVE_C99_INLINE": None}, True),
    "disable-deprecated": ("GSL_DISABLE_DEPRECATED defined in config.h (library only: the tests "
                           "exercise deprecated functions)", [], {"GSL_DISABLE_DEPRECATED": "1"}, False),
    "sanitize":       ("AddressSanitizer + UndefinedBehaviorSanitizer",
                       ["--disable-shared",
                        "CFLAGS=-g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=undefined",
                        "LDFLAGS=-fsanitize=address,undefined"], None, True),
    "clang":          ("CC=clang", ["CC=clang"], None, True),
    "icc":            ("CC=icc -fp-model precise", ["CC=icc", "CFLAGS=-g -O2 -Wall -fp-model precise"], None, True),
}

_I686 = None


def i686_works():
    """can i686-linux-gnu-gcc link a static program, and does it run here?"""
    global _I686
    if _I686 is None:
        import tempfile
        _I686 = False
        if which("i686-linux-gnu-gcc"):
            with tempfile.TemporaryDirectory() as d:
                c, t = os.path.join(d, "t.c"), os.path.join(d, "t")
                with open(c, "w") as f:
                    f.write("#include <math.h>\nint main(void) { return (int) cos(0.0) - 1; }\n")
                _I686 = (subprocess.run(["i686-linux-gnu-gcc", "-static", c, "-o", t, "-lm"],
                                        capture_output=True).returncode == 0 and
                         subprocess.run([t], capture_output=True).returncode == 0)
    return _I686


def build_available(name):
    if name in ("O3-fma", "x86-32") and not IS_X86:
        return "not an x86 machine"
    if name == "x86-32" and not i686_works():
        return "i686-linux-gnu-gcc not found, or a static 32-bit program does not run"
    need = {"c99": "c99", "clang": "clang", "icc": "icc"}
    if name in need and not which(need[name]):
        return "%s not found" % need[name]
    return None


def edit_config_h(path, edits):
    """set (value) or unset (None) macros in a generated config.h"""
    with open(path) as f:
        text = f.read()
    for macro, value in edits.items():
        pat = re.compile(r"^(#define %s\b.*|/\* #undef %s \*/)$" % (macro, macro), re.M)
        new = ("#define %s %s" % (macro, value)) if value is not None else ("/* #undef %s */" % macro)
        if not pat.search(text):
            raise CheckError("config.h has no entry for %s" % macro)
        text = pat.sub(new, text)
    with open(path, "w") as f:
        f.write(text)


def configure_and_build(ctx, name, args, edits=None, env=None, host=None, target_dir=None, src=None):
    """configure and build out of tree (of the revision being checked,
    or of the source tree src); returns the build directory"""
    src = src or prepare_source(ctx)
    bdir = target_dir or os.path.join(ctx.workdir, "build", name)
    if os.path.exists(bdir):
        shutil.rmtree(bdir)
    os.makedirs(bdir)
    cfg = [os.path.join(src, "configure")] + list(args)
    if host:
        cfg.append("--host=%s" % host)
    ctx.run(cfg, bdir, env=env)
    if edits:
        edit_config_h(os.path.join(bdir, "config.h"), edits)
    ctx.make(bdir, env=env)
    return bdir


def built_from(bdir):
    """the source tree a build directory was last built from (its marker), or None"""
    try:
        with open(os.path.join(bdir, ".prerelease-built")) as f:
            return f.read().strip() or None
    except OSError:
        return None


def mark_built(bdir, src):
    with open(os.path.join(bdir, ".prerelease-built"), "w") as f:
        f.write(src + "\n")


def ensure_default_build(ctx):
    """the default build, which several checks use"""
    if "default" in ctx.builds and os.path.exists(ctx.builds["default"]):
        return ctx.builds["default"]
    bdir = os.path.join(ctx.workdir, "build", "default")
    # reuse a default build only if it was built from this revision's export
    if os.path.exists(os.path.join(bdir, ".libs", "libgsl.a")) and built_from(bdir) == prepare_source(ctx):
        ctx.builds["default"] = bdir
        return bdir
    ctx.status("  (building the default configuration first)")
    bdir = configure_and_build(ctx, "default", [])
    mark_built(bdir, prepare_source(ctx))
    ctx.builds["default"] = bdir
    return bdir


WARN_RE = re.compile(r"^(?P<file>[^\s:]+\.[ch]):(?P<line>\d+):(?:\d+:)? warning: (?P<msg>.*)$")


def compiler_warnings(log):
    """unique compiler warnings in a log (file:line: message)"""
    seen = []
    try:
        with open(log, errors="replace") as f:
            for l in f:
                m = WARN_RE.match(l.strip())
                if m:
                    w = "%s:%s: %s" % (os.path.relpath(m.group("file"), "/") if m.group("file").startswith("/")
                                       else m.group("file"), m.group("line"), m.group("msg"))
                    w = re.sub(r"^.*?src-[0-9a-f]+/", "", w)
                    if w not in seen:
                        seen.append(w)
    except OSError:
        pass
    return seen


def test_results(bdir):
    """results of the automake test harness: (n run, failures [(test, lines)])"""
    nrun, fails = 0, []
    for trs in sorted(glob.glob(os.path.join(bdir, "**", "*.trs"), recursive=True)):
        with open(trs, errors="replace") as f:
            res = re.findall(r"^:test-result: (\S+)", f.read(), re.M)
        name = os.path.relpath(trs, bdir)[:-4]
        nrun += 1
        if any(r in ("FAIL", "ERROR", "XPASS") for r in res):
            lines = []
            logf = trs[:-4] + ".log"
            try:
                with open(logf, errors="replace") as f:
                    lines = [l.rstrip() for l in f if l.startswith("FAIL")][:8]
            except OSError:
                pass
            if not lines:
                lines = tail(logf, 5)
            fails.append((name, lines))
    return nrun, fails



def make_var(bdir, var):
    """the value of a variable in a generated Makefile"""
    p = subprocess.run(["make", "-s", "--no-print-directory", "-f", "Makefile", "-f", "-", "prerelease-print"],
                       cwd=bdir, input="prerelease-print:\n\t@echo $(%s)\n" % var,
                       capture_output=True, text=True)
    return p.stdout.split() if p.returncode == 0 else []


def expected_tests(bdir):
    """{directory: [test names]} from the TESTS variables of the generated
    Makefiles ("." is the top-level directory)"""
    tests = {}
    for d in ["."] + make_var(bdir, "SUBDIRS"):
        if d != "." and not os.path.exists(os.path.join(bdir, d, "Makefile")):
            continue
        t = make_var(os.path.join(bdir, d), "TESTS")
        if t:
            tests[d] = t
    return tests


def compiled_tests(tests):
    """the expected tests which are compiled programs (not scripts)"""
    return [(d, t) for d in tests for t in tests[d] if not t.endswith((".sh", ".test"))]


def build_errors(log, n=10):
    """compiler and linker errors in a log"""
    rx = re.compile(r"error:|undefined reference|cannot find -l")
    out = []
    for l in tail(log, 100000):
        if rx.search(l) and "-Werror" not in l and "Error " not in l:
            l = re.sub(r"^.*?src-[0-9a-f]+/", "", l.strip())
            if l not in out:
                out.append(l)
    return out[:n]


def build_test_programs(ctx, bdir, env=None):
    """build the library and every test program without running them;
    returns ({dir: [tests]}, details of build failures)"""
    tests = expected_tests(bdir)
    if not tests:
        raise CheckError("no TESTS found in the Makefiles of %s" % bdir)
    rc = ctx.make(bdir, ["check"], extra=["TESTS="], env=env, check=False, keep_going=True)
    details = []
    if rc != 0:
        details = ["building the test programs failed (exit %d), e.g.:" % rc] + \
                  ["    " + l for l in build_errors(ctx.current_log)]
    return tests, details


# the slowest test directories go first, so that they overlap the others
SLOW_TESTS = ["specfunc", "linalg", "multifit_nlinear", "multilarge_nlinear", "eigen", "randist",
              "rng", "integration", "bspline", "monte", "statistics", "movstat", "ode-initval2", "multifit"]


def run_pool(ctx, items, start, finish):
    """run start(item) -> Popen for each item, ctx.jobs at a time, and call
    finish(item, exit code, seconds) as each one ends"""
    queue, procs = list(items), []
    while queue or procs:
        while queue and len(procs) < ctx.jobs:
            item = queue.pop(0)
            procs.append((start(item), item, time.time()))
        time.sleep(1)
        for entry in list(procs):
            p, item, t0 = entry
            rc = p.poll()
            if rc is None and time.time() - t0 > ctx.args.timeout * 3600:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                rc = -9999
            if rc is not None:
                procs.remove(entry)
                finish(item, rc, time.time() - t0)


def run_tests(ctx, bdir, env=None, prefix="", limit_kb=None):
    """build all the test programs, then run make check in the test
    directories ctx.jobs at a time (each with -j1); returns (status,
    summary, details), counting a test without a result as failed;
    limit_kb limits the virtual memory of the test runs (ulimit -v)"""
    tests, details = build_test_programs(ctx, bdir, env)
    for trs in glob.glob(os.path.join(bdir, "**", "*.trs"), recursive=True):
        os.remove(trs)
    dirs = sorted(tests, key=lambda d: (d not in SLOW_TESTS, SLOW_TESTS.index(d) if d in SLOW_TESTS else 0, d))
    logdir = os.path.join(bdir, "prerelease-check")
    os.makedirs(logdir, exist_ok=True)
    e = dict(os.environ, LC_ALL="C")
    if env:
        e.update(env)
    ctx.status("  %srunning the tests of %d directories (%d at a time)" % (prefix, len(dirs), ctx.jobs))

    def start(d):
        # at the top level, check-TESTS does not recurse into the subdirectories
        target = "check-TESTS" if d == "." else "check"
        out = open(os.path.join(logdir, d.replace("/", "_") + ".log"), "w")
        def preexec():
            ctx._preexec()
            if limit_kb:
                resource.setrlimit(resource.RLIMIT_AS, (limit_kb * 1024, limit_kb * 1024))
        return subprocess.Popen(["make", "-j1", "-k", target], cwd=os.path.join(bdir, d), env=e,
                                stdout=out, stderr=subprocess.STDOUT, preexec_fn=preexec)

    def finish(d, rc, secs):
        if rc == -9999:
            details.append("%s: make check timed out" % d)
        if secs > 300:
            ctx.status("    %s: %s" % (d, human(secs)))

    run_pool(ctx, dirs, start, finish)
    nrun, fails = test_results(bdir)
    # automake names the results without the .test extension (pkgconfig.test -> pkgconfig.trs)
    missing = [os.path.normpath(os.path.join(d, t)) for d in tests for t in tests[d]
               if not os.path.exists(os.path.join(bdir, d, re.sub(r"\.test$", "", t) + ".trs"))]
    for t, lines in fails:
        details.append("FAIL %s" % t)
        details += ["    " + l for l in lines]
    if missing:
        details.append("%d tests have no result (did not build or did not run): %s"
                       % (len(missing), " ".join(missing)))
    nexp = sum(len(v) for v in tests.values())
    status = "FAIL" if fails or missing or details else "PASS"
    summary = "%d/%d tests run, %d failed" % (nrun, nexp, len(fails))
    if missing:
        summary += ", %d MISSING" % len(missing)
    return status, summary, details

def check_build(ctx, name):
    desc, args, edits, run_check = BUILDS[name]
    why = build_available(name)
    if why:
        return Result("SKIP", why)
    bdir = configure_and_build(ctx, name, args, edits)
    details = []
    warns = compiler_warnings(ctx.current_log)
    if run_check:
        env = None
        if name == "sanitize":
            env = {"ASAN_OPTIONS": "detect_leaks=0:abort_on_error=0", "UBSAN_OPTIONS": "print_stacktrace=1"}
        status, summary, details = run_tests(ctx, bdir, env)
        if name == "sanitize":
            san = grep_logs(bdir, r"runtime error:|ERROR: AddressSanitizer")[:10]
            if san:
                details += ["sanitizer reports:"] + ["    " + l for l in san]
                status = "FAIL"
    else:
        status, summary = "PASS", "library built"
    implicit = [w for w in warns if "implicit declaration" in w]
    if implicit:
        # an undeclared function is assumed to return int: miscompiled code
        details = ["%d implicit declarations (the functions are assumed to return int):" % len(implicit)] + \
                  ["    " + w for w in implicit[:20]] + details
        status = "FAIL"
        summary += ", %d IMPLICIT DECLARATIONS" % len(implicit)
    if warns:
        details += ["%d compiler warnings, e.g.:" % len(warns)] + ["    " + w for w in warns[:10]]
        status = worst(status, "WARN")
        summary += ", %d compiler warnings" % len(warns)
    if name == "default":
        mark_built(bdir, prepare_source(ctx))
        ctx.builds["default"] = bdir
    elif status == "PASS" and not ctx.args.keep_builds:
        shutil.rmtree(bdir, ignore_errors=True)
    return Result(status, summary, details)


def grep_logs(bdir, pattern):
    rx = re.compile(pattern)
    out = []
    for logf in sorted(glob.glob(os.path.join(bdir, "**", "*.log"), recursive=True)):
        try:
            with open(logf, errors="replace") as f:
                for l in f:
                    if rx.search(l):
                        out.append("%s: %s" % (os.path.relpath(logf, bdir), l.strip()))
        except OSError:
            pass
    return out


# ---------------------------------------------------------------------------
# static checks


def check_symbols(ctx):
    """every exported symbol starts with gsl_/GSL_ (cblas_ in libgslcblas)"""
    bdir = ensure_default_build(ctx)
    bad = []
    for lib, prefixes in ((".libs/libgsl.a", ("gsl_", "GSL_")),
                          ("cblas/.libs/libgslcblas.a", ("cblas_", "gsl_", "GSL_"))):
        path = os.path.join(bdir, lib)
        out = subprocess.run(["nm", "-A", "-g", "-P", "--defined-only", path],
                             capture_output=True, text=True).stdout
        for line in out.splitlines():
            parts = line.split()
            if len(parts) < 3:
                continue
            sym, typ = parts[1], parts[2]
            if typ in ("U", "w", "v"):
                continue
            # as in HACKING, the prefix may follow underscores (_gsl_sf_..._cs tables)
            if not sym.lstrip("_").startswith(prefixes):
                obj = parts[0].split("[")[-1].rstrip(":]")
                bad.append("%s: %s (%s) in %s" % (lib, sym, typ, obj))
    if bad:
        return Result("FAIL", "%d exported symbols without the gsl_ prefix" % len(bad), bad[:40])
    return Result("PASS", "all exported symbols prefixed")


def check_config_h(ctx):
    """installed headers must not use config.h; compiled library sources must"""
    src = prepare_source(ctx)
    headers = sorted(glob.glob(os.path.join(src, "*", "gsl_*.h")) + glob.glob(os.path.join(src, "gsl_*.h")))
    bad_h = [os.path.relpath(h, src) for h in headers if "config.h" in open(h, errors="replace").read()]
    # the sources compiled into the libraries (lib*_la_SOURCES); other .c
    # files are #included by these, and see config.h through them
    compiled = set()
    for mk in glob.glob(os.path.join(src, "*", "Makefile.am")):
        d = os.path.relpath(os.path.dirname(mk), src)
        text = open(mk, errors="replace").read().replace("\\\n", " ")
        for m in re.finditer(r"^lib\w+_la_SOURCES\s*=(.*)$", text, re.M):
            compiled.update(os.path.join(d, f) for f in m.group(1).split() if f.endswith(".c"))
    missing = [c for c in sorted(compiled) if os.path.exists(os.path.join(src, c))
               and "config.h" not in open(os.path.join(src, c), errors="replace").read()]
    details = ["installed header uses config.h: %s" % h for h in bad_h]
    cblas = [m for m in missing if m.startswith("cblas/")]
    other = [m for m in missing if not m.startswith("cblas/")]
    if other:
        details.append("%d library sources do not include config.h:" % len(other))
        details += ["    " + m for m in other]
    if cblas:
        details.append("%d cblas sources do not include config.h" % len(cblas))
    status = "FAIL" if bad_h else ("WARN" if missing else "PASS")
    return Result(status, "%d headers use config.h, %d of %d library sources without config.h"
                  % (len(bad_h), len(missing), len(compiled)), details)


def check_headers(ctx):
    """all headers together (C, C++), and each header on its own in several modes"""
    bdir = ensure_default_build(ctx)
    gsldir = os.path.join(bdir, "gsl")
    headers = sorted(os.path.basename(h) for h in glob.glob(os.path.join(gsldir, "*.h")))
    tdir = os.path.join(ctx.workdir, "headers")
    os.makedirs(tdir, exist_ok=True)
    inc = ["-I", bdir]
    failures = []

    # all headers concatenated (HACKING)
    allc = os.path.join(tdir, "all.c")
    with open(allc, "w") as f:
        for h in headers:
            f.write(open(os.path.join(gsldir, h), errors="replace").read() + "\n")
    for comp, flags in (("gcc", ["-g", "-Wall", "-W"]), ("g++", ["-g", "-Wall", "-W", "-x", "c++"])):
        rc = ctx.run([comp] + flags + inc + ["-c", allc, "-o", os.devnull], tdir, check=False)
        if rc != 0:
            failures.append("all headers concatenated with %s" % comp)

    # each header on its own
    modes = [
        ("c89", ["gcc", "-std=c89", "-pedantic", "-Wall", "-W", "-Werror"], ""),
        ("c99+inline", ["gcc", "-std=c99", "-pedantic", "-Wall", "-W", "-Werror", "-DHAVE_INLINE"], ""),
        ("c11+complex.h", ["gcc", "-std=c11", "-pedantic", "-Wall", "-W", "-Werror"], "#include <complex.h>\n"),
        ("c89+disable-deprecated", ["gcc", "-std=c89", "-pedantic", "-Wall", "-Werror", "-DGSL_DISABLE_DEPRECATED"], ""),
        ("c++", ["g++", "-x", "c++", "-Wall", "-W", "-Werror"], ""),
    ]
    for h in headers:
        for mode, cmd, pre in modes:
            t = os.path.join(tdir, "one.c")
            with open(t, "w") as f:
                f.write("%s#include <gsl/%s>\nint gsl_header_check_dummy;\n" % (pre, h))
            rc = ctx.run(cmd + inc + ["-c", t, "-o", os.devnull], tdir, check=False)
            if rc != 0:
                failures.append("%s alone in %s" % (h, mode))
    if failures:
        return Result("FAIL", "%d header compilations failed" % len(failures),
                      failures[:40] + ["see %s" % ctx.current_log])
    return Result("PASS", "%d headers, together and alone in %d modes" % (len(headers), len(modes)))


def check_version(ctx):
    """version strings consistent across configure.ac, NEWS and doc/conf.py"""
    src = prepare_source(ctx)
    cfg = open(os.path.join(src, "configure.ac")).read()
    details, status = [], "PASS"
    m = re.search(r"AC_INIT\(\[gsl\],\s*\[([^\]]+)\]", cfg)
    version = m.group(1) if m else "?"
    if "+" in version:
        details.append("AC_INIT version is %s: remove the '+' for a release" % version)
        status = "WARN"
    v = version.rstrip("+")
    news = open(os.path.join(src, "NEWS")).readline().strip()
    if not news.startswith("* What is new in gsl-%s:" % v):
        details.append("NEWS starts with %r, expected '* What is new in gsl-%s:'" % (news, v))
        status = "WARN"
    conf = open(os.path.join(src, "doc", "conf.py")).read()
    for key in ("release", "version"):
        mm = re.search(r"^%s\s*=\s*'([^']*)'" % key, conf, re.M)
        if mm and mm.group(1) != v:
            details.append("doc/conf.py %s = '%s', expected '%s'" % (key, mm.group(1), v))
            status = "WARN"
    year = time.strftime("%Y")
    mm = re.search(r"^copyright\s*=\s*'([^']*)'", conf, re.M)
    if mm and year not in mm.group(1):
        details.append("doc/conf.py copyright '%s' does not include %s" % (mm.group(1), year))
        status = "WARN"
    lt = {k: re.search(r"^%s=(\d+)" % k, cfg, re.M) for k in ("GSL_CURRENT", "GSL_REVISION", "GSL_AGE")}
    lts = ":".join(lt[k].group(1) if lt[k] else "?" for k in ("GSL_CURRENT", "GSL_REVISION", "GSL_AGE"))
    return Result(status, "version %s, libtool %s" % (version, lts), details)


def check_autoreconf(ctx):
    """autoreconf -W all: obsolete constructs in configure.ac and Makefile.am"""
    src = prepare_source(ctx)
    copy = os.path.join(ctx.workdir, "src-autoreconf")
    if os.path.exists(copy):
        shutil.rmtree(copy)
    shutil.copytree(src, copy, symlinks=True)
    ctx.run(["autoreconf", "-W", "all", "-f"], copy, check=False)
    warns = [l for l in tail(ctx.current_log, 2000) if "warning" in l.lower()]
    shutil.rmtree(copy, ignore_errors=True)
    if warns:
        return Result("WARN", "%d autoreconf warnings" % len(warns), warns[:30])
    return Result("PASS", "no warnings")


# ---------------------------------------------------------------------------
# dynamic checks


def check_verbose(ctx):
    """make check with GSL_TEST_VERBOSE=1 (format string errors in verbose output)"""
    bdir = ensure_default_build(ctx)
    status, summary, details = run_tests(ctx, bdir, env={"GSL_TEST_VERBOSE": "1"})
    return Result(status, summary, details)


def check_examples(ctx):
    """build the documentation examples and compare their output with the references"""
    bdir = ensure_default_build(ctx)
    src = prepare_source(ctx)
    ctx.make(os.path.join(bdir, "doc", "examples"), ["check"], check=False, keep_going=True)
    rc = ctx.run([sys.executable, os.path.join(src, "scripts", "test_examples.py"),
                  "--exdir", os.path.join(bdir, "doc", "examples"),
                  "--refdir", os.path.join(src, "doc", "examples")], src, check=False)
    lines = [l for l in tail(ctx.current_log, 400) if l.startswith(("FAIL", "ERROR", "MISSING"))]
    summ = [l for l in tail(ctx.current_log, 10) if re.search(r"\d+ (passed|failed)|PASS|FAIL", l)]
    return Result("FAIL" if rc != 0 else "PASS", summ[-1] if summ else "exit %d" % rc, lines[:30])


def check_valgrind(ctx):
    """every test program under valgrind memcheck (static build, -g, no optimization)"""
    if not which("valgrind"):
        return Result("SKIP", "valgrind not found")
    # long double padding written by the fwrite functions (see the file)
    supp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "valgrind.supp")
    # static, so that the test programs are executables and not libtool wrapper scripts
    bdir = configure_and_build(ctx, "valgrind", ["--disable-shared", "CFLAGS=-g -Wall"])
    return run_programs(ctx, bdir, "valgrind",
                        lambda prog: ["valgrind", "--leak-check=full", "--error-exitcode=100",
                                      "--errors-for-leak-kinds=definite", "--suppressions=" + supp,
                                      "./" + os.path.basename(prog)],
                        summary_re=r"ERROR SUMMARY|definitely lost")


def run_programs(ctx, bdir, label, command, env=None, summary_re=None):
    """build the test programs of a static build and run each one as
    command(program), ctx.jobs at a time; a test which did not build fails"""
    tests, details = build_test_programs(ctx, bdir, env)
    expected = [os.path.join(bdir, d, t) for d, t in compiled_tests(tests)]
    progs = [p for p in expected if is_elf(p)]
    missing = [os.path.relpath(p, bdir) for p in expected if not is_elf(p)]
    odir = os.path.join(bdir, "prerelease-" + label)
    os.makedirs(odir, exist_ok=True)
    e = dict(os.environ, LC_ALL="C")
    if env:
        e.update(env)
    ctx.status("  running %d test programs under %s (%d at a time)" % (len(progs), label, ctx.jobs))
    results = []

    def out_file(prog):
        return os.path.join(odir, os.path.relpath(prog, bdir).replace("/", "_") + ".txt")

    def start(prog):
        return subprocess.Popen(command(prog), cwd=os.path.dirname(prog), env=e, stdout=open(out_file(prog), "w"),
                                stderr=subprocess.STDOUT, preexec_fn=ctx._preexec)

    def finish(prog, rc, secs):
        out = out_file(prog)
        with open(out, errors="replace") as f:
            text = f.readlines()
        fl = [l.rstrip() for l in text if l.startswith("FAIL")]
        extra = [l.strip() for l in text if summary_re and re.search(summary_re, l)][-2:]
        name = os.path.relpath(prog, bdir)
        results.append((name, rc, fl, extra, out))
        ctx.status("    %s %d/%d %s: %s" % (label, len(results), len(progs), name,
                                            "ok" if rc == 0 and not fl else "exit %d, %d FAIL" % (rc, len(fl))))

    run_pool(ctx, progs, start, finish)
    bad = [r for r in results if r[1] != 0 or r[2]]
    for name, rc, fl, extra, out in bad:
        sig = " (signal %d)" % (rc - 128) if rc > 128 else (" (killed by signal %d)" % -rc if rc < 0 else "")
        details.append("%s: exit %d%s, %d FAIL lines, see %s" % (name, rc, sig, len(fl), out))
        details += ["    " + re.sub(r" \[\d+\]$", "", l) for l in fl[:5] + extra]
    if missing:
        details.append("%d test programs did not build: %s" % (len(missing), " ".join(missing)))
    status = "FAIL" if bad or missing or details else "PASS"
    summary = "%d/%d test programs run, %d failed" % (len(results), len(expected), len(bad))
    if missing:
        summary += ", %d MISSING" % len(missing)
    if status == "PASS" and not ctx.args.keep_builds:
        shutil.rmtree(bdir, ignore_errors=True)
    return Result(status, summary, details)


# ---------------------------------------------------------------------------
# distribution


def check_distcheck(ctx, install_sh=False):
    bdir = ensure_default_build(ctx)
    extra = []
    if install_sh:
        extra = ["INSTALL=%s -c" % os.path.join(prepare_source(ctx), "install-sh")]
    rc = ctx.make(bdir, ["distcheck"], extra=extra, check=False)
    ok = any("archives ready for distribution" in l for l in tail(ctx.current_log, 30))
    if rc != 0 or not ok:
        return Result("FAIL", "make distcheck failed", ["see %s" % ctx.current_log] + tail(ctx.current_log, 20))
    return Result("PASS", "distribution tarball builds and passes make check")


def check_tarball(ctx):
    """make dist, unpack, configure, then build the library and the test
    programs under ulimit -v (HACKING: compiling must not need excessive
    memory), and run the tests under a looser limit"""
    bdir = ensure_default_build(ctx)
    ctx.make(bdir, ["dist-gzip"])
    tars = sorted(glob.glob(os.path.join(bdir, "gsl-*.tar.gz")), key=os.path.getmtime)
    if not tars:
        raise CheckError("make dist produced no tarball")
    tdir = os.path.join(ctx.workdir, "tarball")
    if os.path.exists(tdir):
        shutil.rmtree(tdir)
    os.makedirs(tdir)
    ctx.run(["tar", "-xzf", tars[-1]], tdir)
    top = os.path.join(tdir, os.path.basename(tars[-1])[:-len(".tar.gz")])
    ctx.run(["./configure"], top)
    kb, test_kb = ctx.args.ulimit_kb, ctx.args.test_ulimit_kb
    name = os.path.basename(tars[-1])
    for what, targets in (("the library", ""), ("the test programs", " -k check TESTS=")):
        rc = ctx.run("ulimit -v %d && make -j%d%s" % (kb, ctx.jobs, targets), top, shell=True, check=False)
        if rc != 0:
            return Result("FAIL", "%s: building %s failed under ulimit -v %d" % (name, what, kb),
                          ["see %s" % ctx.current_log] + build_errors(ctx.current_log) +
                          [l for l in tail(ctx.current_log, 400) if "memory" in l.lower()][:10])
    status, summary, details = run_tests(ctx, top, limit_kb=test_kb)
    if not ctx.args.keep_builds and status == "PASS":
        shutil.rmtree(tdir, ignore_errors=True)
    return Result(status, "%s built under ulimit -v %d, tests under %d: %s" % (name, kb, test_kb, summary), details)


# ---------------------------------------------------------------------------
# documentation


def check_docs(ctx, target):
    if not which("sphinx-build"):
        return Result("SKIP", "sphinx-build not found")
    if target == "latexpdf" and not which("latexmk"):
        return Result("SKIP", "latexmk not found")
    if target == "info" and not which("makeinfo"):
        return Result("SKIP", "makeinfo not found")
    src = prepare_source(ctx)
    out = os.path.join(ctx.workdir, "docs")
    rc = ctx.run(["sphinx-build", "-M", target, os.path.join(src, "doc"), out, "-q"], src, check=False)
    log = tail(ctx.current_log, 5000)
    warns = [l for l in log if "WARNING" in l or "ERROR" in l]
    latex_err = [l for l in log if l.startswith("!")]
    details = warns[:25] + latex_err[:10]
    if rc != 0:
        return Result("FAIL", "sphinx-build -M %s failed (exit %d)" % (target, rc), details + tail(ctx.current_log, 10))
    return Result("WARN" if warns else "PASS", "%d warnings" % len(warns), details)


# ---------------------------------------------------------------------------
# ABI


def previous_release(ctx):
    if ctx.args.abi_base:
        return ctx.args.abi_base
    tags = subprocess.run(["git", "tag", "--list", "release-*", "--merged", ctx.args.rev], cwd=ctx.repo,
                          capture_output=True, text=True).stdout.split()

    def key(t):
        return [int(x) for x in re.findall(r"\d+", t)]
    tags = sorted(tags, key=key)
    return tags[-1] if tags else None


def lt_version(configure_ac):
    text = open(configure_ac).read()
    return tuple(int(re.search(r"^%s=(\d+)" % k, text, re.M).group(1))
                 for k in ("GSL_CURRENT", "GSL_REVISION", "GSL_AGE"))


def exported(lib):
    out = subprocess.run(["nm", "-D", "--defined-only", lib], capture_output=True, text=True).stdout
    return set(l.split()[-1] for l in out.splitlines() if len(l.split()) >= 3 and l.split()[1] in "TDRBVW")


def check_abi(ctx):
    """interfaces removed/added/changed since the previous release, and the libtool numbers"""
    base = previous_release(ctx)
    if not base:
        return Result("SKIP", "no previous release tag found (use --abi-base)")
    src = prepare_source(ctx)
    bdir = ensure_default_build(ctx)
    # the previous release, built shared with debug information
    bsrc = os.path.join(ctx.workdir, "abi-base-src")
    bbld = os.path.join(ctx.workdir, "abi-base-build")
    # reuse the base build only if it is of the same release
    if built_from(bbld) != base:
        ctx.status("  building %s for comparison" % base)
        for d in (bsrc, bbld):
            if os.path.exists(d):
                shutil.rmtree(d)
        os.makedirs(bsrc)
        ctx.run("git archive %s | tar -x -C %s" % (shlex.quote(base), shlex.quote(bsrc)), ctx.repo, shell=True)
        ctx.run(["sh", "autogen.sh"], bsrc)
        configure_and_build(ctx, "abi-base", ["CFLAGS=-g -O2"], target_dir=bbld, src=bsrc)
        mark_built(bbld, base)

    details, status = [], "PASS"
    removed_total, added_total = 0, 0
    for lib in (".libs/libgsl.so", "cblas/.libs/libgslcblas.so"):
        old, new = os.path.join(bbld, lib), os.path.join(bdir, lib)
        if not (os.path.exists(old) and os.path.exists(new)):
            return Result("FAIL", "shared library %s missing (configure with shared libraries)" % lib)
        so, sn = exported(old), exported(new)
        removed, added = sorted(so - sn), sorted(sn - so)
        removed_total += len(removed) if "cblas" not in lib else 0
        added_total += len(added) if "cblas" not in lib else 0
        details.append("%s: %d exported symbols removed, %d added" % (os.path.basename(lib), len(removed), len(added)))
        details += ["    removed: %s" % s for s in removed[:30]]
        if which("abidiff"):
            out = os.path.join(ctx.workdir, "abidiff-%s.txt" % os.path.basename(lib))
            with open(out, "w") as f:
                subprocess.run(["abidiff", "--headers-dir1", os.path.join(bbld, "gsl"),
                                "--headers-dir2", os.path.join(bdir, "gsl"), old, new],
                               stdout=f, stderr=subprocess.STDOUT)
            txt = open(out, errors="replace").read()
            changed = re.findall(r"^  \[C\] '([^']*)'", txt, re.M)
            sizes = re.findall(r"underlying type '([^']*)' at [^\n]*\n\s*type size changed from (\d+) to (\d+)", txt)
            details.append("    abidiff: %d changed interfaces (full report: %s)" % (len(changed), out))
            details += ["      changed: %s" % c for c in changed[:15]]
            details += ["      type %s: size %s -> %s bits" % s for s in sizes[:15]]
            if changed:
                status = worst(status, "WARN")
    # libtool numbers
    oc, orv, oa = lt_version(os.path.join(bsrc, "configure.ac"))
    nc, nrv, na = lt_version(os.path.join(src, "configure.ac"))
    details.append("libtool version %s: %d:%d:%d, now: %d:%d:%d" % (base, oc, orv, oa, nc, nrv, na))
    if removed_total:
        expect = "%d:0:0 (interfaces removed)" % (oc + 1)
        okv = (nc, nrv, na) == (oc + 1, 0, 0)
    elif added_total:
        expect = "%d:0:%d (interfaces added; age 0 instead if any change is incompatible)" % (oc + 1, oa + 1)
        okv = (nc, nrv) == (oc + 1, 0) and na in (oa + 1, 0)
    else:
        expect = "%d:%d:%d (no interface change)" % (oc, orv + 1, oa)
        okv = (nc, nrv, na) == (oc, orv + 1, oa)
    details.append("expected %s" % expect)
    if not okv:
        status = "FAIL"
    return Result(status, "vs %s: libgsl %d symbols removed, %d added; libtool %d:%d:%d %s"
                  % (base, removed_total, added_total, nc, nrv, na, "ok" if okv else "INCONSISTENT"), details)


# ---------------------------------------------------------------------------
# QEMU cross-architecture runs

QEMU_TARGETS = {
    "aarch64-linux-gnu": "qemu-aarch64-static",
    "powerpc64le-linux-gnu": "qemu-ppc64le-static",
    "powerpc64-linux-gnu": "qemu-ppc64-static",
    "sparc64-linux-gnu": "qemu-sparc64-static",
    "arm-linux-gnueabihf": "qemu-arm-static",
    "riscv64-linux-gnu": "qemu-riscv64-static",
    "s390x-linux-gnu": "qemu-s390x-static",
    "mips64el-linux-gnuabi64": "qemu-mips64el-static",
    "i686-linux-gnu": "qemu-i386-static",
}


def qemu_available(target):
    q = QEMU_TARGETS[target]
    if not which(target + "-gcc"):
        return "%s-gcc not found" % target
    if not which(q):
        return "%s not found" % q
    return None


def check_qemu(ctx, target):
    why = qemu_available(target)
    if why:
        return Result("SKIP", why)
    q = QEMU_TARGETS[target]
    bdir = configure_and_build(ctx, "qemu-" + target, ["--disable-shared"], host=target)
    return run_programs(ctx, bdir, target, lambda prog: [q, "./" + os.path.basename(prog)],
                        env={"QEMU_LD_PREFIX": "/usr/%s" % target})


# ---------------------------------------------------------------------------
# check registry

def all_checks():
    """(name, group, description, function); in the order they run"""
    checks = [
        ("version", "static", "version numbers in configure.ac, NEWS, doc/conf.py", check_version),
        ("config-h", "static", "config.h use in headers and sources", check_config_h),
        ("build:default", "build", BUILDS["default"][0], lambda c: check_build(c, "default")),
        ("symbols", "static", "exported symbols prefixed with gsl_", check_symbols),
        ("headers", "static", "headers compiled together and alone (C89, C99, C11, C++)", check_headers),
        ("examples", "dynamic", "documentation examples (scripts/test_examples.py)", check_examples),
        ("verbose", "dynamic", "make check with GSL_TEST_VERBOSE=1", check_verbose),
        ("abi", "abi", "interfaces and libtool version vs the previous release", check_abi),
        ("autoreconf", "static", "autoreconf -W all", check_autoreconf),
        ("distcheck", "dist", "make distcheck", lambda c: check_distcheck(c)),
        ("distcheck-install-sh", "dist", "make distcheck INSTALL=install-sh", lambda c: check_distcheck(c, True)),
        ("tarball", "dist", "make dist, unpack, build under ulimit -v, make check", check_tarball),
    ]
    for name in BUILDS:
        if name != "default":
            checks.append(("build:" + name, "build", BUILDS[name][0], (lambda n: lambda c: check_build(c, n))(name)))
    checks += [
        ("docs:html", "docs", "Sphinx HTML manual", lambda c: check_docs(c, "html")),
        ("docs:pdf", "docs", "Sphinx PDF manual (LaTeX)", lambda c: check_docs(c, "latexpdf")),
        ("docs:info", "docs", "Sphinx info manual", lambda c: check_docs(c, "info")),
    ]
    for t in QEMU_TARGETS:
        checks.append(("qemu:" + t, "qemu", "test programs under QEMU (%s)" % t, (lambda tt: lambda c: check_qemu(c, tt))(t)))
    checks.append(("valgrind", "dynamic", "every test program under valgrind memcheck", check_valgrind))
    return checks


QUICK = ["version", "config-h", "build:default", "symbols", "headers", "examples", "abi", "docs:html"]


def select(checks, args):
    names = [c[0] for c in checks]
    sel = names
    if args.quick:
        sel = [n for n in names if n in QUICK]
    if args.only:
        pats = args.only.split(",")
        sel = [n for n in names if any(fnmatch.fnmatch(n, p) or n.split(":")[0] == p for p in pats)]
    if args.skip:
        pats = args.skip.split(",")
        sel = [n for n in sel if not any(fnmatch.fnmatch(n, p) or n.split(":")[0] == p for p in pats)]
    return [c for c in checks if c[0] in sel]


# ---------------------------------------------------------------------------
# summary


def write_summary(ctx, results, durations, started):
    order = {"FAIL": 0, "WARN": 1, "PASS": 2, "SKIP": 3}
    lines = ["GSL pre-release checks: %s (%s)" % (ctx.args.rev, ctx.commit or "?"),
             "work directory: %s" % ctx.workdir,
             "started %s, %s elapsed" % (time.strftime("%Y-%m-%d %H:%M", time.localtime(started)),
                                       human(time.time() - started)),
             ""]
    counts = {}
    for name, r in results.items():
        counts[r.status] = counts.get(r.status, 0) + 1
    lines.append("  ".join("%s %d" % (s, counts.get(s, 0)) for s in ("FAIL", "WARN", "PASS", "SKIP")))
    lines.append("")
    w = max([len(n) for n in results] + [10])
    for name, r in results.items():
        lines.append("%-4s  %-*s  %8s  %s" % (r.status, w, name, human(durations.get(name, 0)), r.summary))
    for name, r in sorted(results.items(), key=lambda kv: order[kv[1].status]):
        if r.status in ("FAIL", "WARN") and r.details:
            lines += ["", "=== %s: %s" % (r.status, name)] + ["  " + d for d in r.details]
    lines += ["", "per-check logs: %s/logs/" % ctx.workdir]
    text = "\n".join(lines) + "\n"
    with open(os.path.join(ctx.workdir, "summary.txt"), "w") as f:
        f.write(text)
    with open(os.path.join(ctx.workdir, "summary.json"), "w") as f:
        json.dump({"rev": ctx.args.rev, "commit": ctx.commit,
                   "results": {n: dict(r.to_dict(), seconds=durations.get(n, 0)) for n, r in results.items()}},
                  f, indent=1)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--workdir", default="prerelease", help="work directory (default ./prerelease)")
    ap.add_argument("--rev", default="HEAD", help="git revision to check (default HEAD)")
    ap.add_argument("--jobs", "-j", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="parallel make jobs and parallel test programs (default: half the CPUs)")
    ap.add_argument("--nice", type=int, default=0, help="run at this nice level (e.g. 19)")
    ap.add_argument("--only", help="comma-separated checks or groups to run (wildcards allowed)")
    ap.add_argument("--skip", help="comma-separated checks or groups to skip (wildcards allowed)")
    ap.add_argument("--quick", action="store_true", help="only the quick checks: %s" % ", ".join(QUICK))
    ap.add_argument("--resume", action="store_true", help="skip checks which passed in a previous run")
    ap.add_argument("--keep-builds", action="store_true", help="keep the build trees of passed checks")
    ap.add_argument("--abi-base", help="release to compare the ABI with (default: latest release-* tag)")
    # HACKING suggests 65536, but GCC 13 needs more than 128 MB of address
    # space for the largest sources (integration/lebedev.c, specfunc/test_sf.c)
    ap.add_argument("--ulimit-kb", type=int, default=262144,
                    help="virtual memory limit for building the tarball (default 262144)")
    ap.add_argument("--test-ulimit-kb", type=int, default=524288,
                    help="virtual memory limit for the tarball make check (default 524288)")
    ap.add_argument("--timeout", type=float, default=8.0, help="time limit per command in hours (default 8)")
    ap.add_argument("--list", action="store_true", help="list the checks and exit")
    args = ap.parse_args()

    checks = all_checks()
    if args.list:
        for name, group, desc, _ in checks:
            why = None
            if name.startswith("build:"):
                why = build_available(name[6:])
            elif name.startswith("qemu:"):
                why = qemu_available(name[5:])
            print("%-26s %-8s %s%s" % (name, group, desc, ("  [unavailable: %s]" % why) if why else ""))
        return 0

    if args.nice:
        os.nice(args.nice)
    ctx = Context(args)
    logdir = os.path.join(ctx.workdir, "logs")
    os.makedirs(logdir, exist_ok=True)
    statefile = os.path.join(ctx.workdir, "state.json")
    state = {}
    if args.resume and os.path.exists(statefile):
        state = json.load(open(statefile))

    selected = select(checks, args)
    results, durations = {}, {}
    started = time.time()
    ctx.status("GSL pre-release checks, %d checks selected, work directory %s" % (len(selected), ctx.workdir))
    try:
        prepare_source(ctx)
        if state.get("commit") and state["commit"] != ctx.commit:
            ctx.status("note: --resume state is for %s, this is %s: running everything" % (state["commit"], ctx.commit))
            state = {}
        state["commit"] = ctx.commit
        state.setdefault("results", {})
        for i, (name, group, desc, fn) in enumerate(selected, 1):
            prev = state["results"].get(name)
            if args.resume and prev and prev["status"] in ("PASS", "SKIP"):
                results[name] = Result.from_dict(prev)
                durations[name] = prev.get("seconds", 0)
                ctx.status("[%d/%d] %s: %s (previous run)" % (i, len(selected), name, prev["status"]))
                continue
            ctx.current_log = os.path.join(logdir, name.replace(":", "_") + ".log")
            with open(ctx.current_log, "w") as f:
                f.write("# %s: %s\n" % (name, desc))
            ctx.status("[%d/%d] %s: %s ..." % (i, len(selected), name, desc))
            t0 = time.time()
            try:
                r = fn(ctx)
            except CheckError as e:
                r = Result("FAIL", str(e), e.details)
            except Exception as e:     # a bug in this script should not stop the run
                r = Result("FAIL", "internal error: %r" % e)
            durations[name] = time.time() - t0
            results[name] = r
            ctx.status("[%d/%d] %s: %s (%s) %s" % (i, len(selected), name, r.status, human(durations[name]), r.summary))
            state["results"][name] = dict(r.to_dict(), seconds=durations[name])
            json.dump(state, open(statefile, "w"), indent=1)
            write_summary(ctx, results, durations, started)
    except KeyboardInterrupt:
        ctx.status("interrupted")
    except CheckError as e:
        ctx.status("cannot prepare the source tree: %s" % e)
        for d in e.details:
            ctx.status("  " + d)
        return 2
    text = write_summary(ctx, results, durations, started)
    print("\n" + text)
    return 1 if any(r.status == "FAIL" for r in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
