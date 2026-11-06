#!/usr/bin/env python3
"""
Run the example programs of the GSL documentation and compare their
output with the reference output files in doc/examples.

The purpose is to detect any change in the output of the examples caused
by a change in the library. Most examples are compared exactly. Examples
which solve ill-conditioned problems, whose last digits depend on the
compiler, platform and floating-point contraction (FMA), are compared
numerically with a tolerance: every number in the output must agree with
the reference within

    |actual - expected| <= atol + rtol * max(|actual|, |expected|)

and all text between the numbers must agree exactly (up to whitespace).

The examples must be built first:

    make -C doc/examples check

Usage:

    scripts/test_examples.py              run all tests
    scripts/test_examples.py -k nlfit     run tests whose name contains "nlfit"
    scripts/test_examples.py --measure    report the largest numeric difference
                                          from the reference for every test
    scripts/test_examples.py --update -k NAME
                                          overwrite the reference files of the
                                          selected tests with the current output
"""

import argparse
import difflib
import fnmatch
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# comparison of outputs
# ---------------------------------------------------------------------------

# numbers: integers, decimals, exponent notation, inf and nan; a number must
# not be part of an identifier such as x1 or a longer number
NUMBER_RE = re.compile(
    r"(?<![A-Za-z_\d.])"
    r"([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
    r"|[-+]?(?:inf(?:inity)?|nan)(?![A-Za-z_]))",
    re.IGNORECASE,
)


def split_numbers(line):
    """Split a line into its text skeleton (whitespace normalized) and numbers."""
    texts = []
    numbers = []
    pos = 0
    for m in NUMBER_RE.finditer(line):
        texts.append(line[pos:m.start()])
        numbers.append(float(m.group(1)))
        pos = m.end()
    texts.append(line[pos:])
    skeleton = "\x00".join(" ".join(t.split()) for t in texts)
    return skeleton, numbers


def number_error(a, b):
    """Return (absolute error, relative error) between two numbers."""
    if a == b or (math.isnan(a) and math.isnan(b)):
        return 0.0, 0.0
    if math.isinf(a) or math.isinf(b) or math.isnan(a) or math.isnan(b):
        return math.inf, math.inf
    abserr = abs(a - b)
    return abserr, abserr / max(abs(a), abs(b))


@dataclass(frozen=True)
class Exact:
    """Outputs must be identical."""

    def describe(self):
        return "exact"

    def compare(self, expected, actual):
        if expected == actual:
            return []
        diff = difflib.unified_diff(expected.splitlines(), actual.splitlines(),
                                    "expected", "actual", lineterm="", n=0)
        return list(diff)[2:]


@dataclass(frozen=True)
class Numeric:
    """Text must be identical, numbers must agree within atol + rtol * max(|a|,|b|)."""

    rtol: float
    atol: float = 0.0

    def describe(self):
        return f"numeric rtol={self.rtol:g} atol={self.atol:g}"

    def compare(self, expected, actual):
        elines = expected.splitlines()
        alines = actual.splitlines()
        if len(elines) != len(alines):
            return [f"number of lines differs: expected {len(elines)}, actual {len(alines)}"]

        errors = []
        for lineno, (eline, aline) in enumerate(zip(elines, alines), start=1):
            eskel, enums = split_numbers(eline)
            askel, anums = split_numbers(aline)
            if eskel != askel or len(enums) != len(anums):
                errors.append(f"line {lineno}: text differs\n  expected: {eline}\n  actual:   {aline}")
                continue
            for a, b in zip(enums, anums):
                abserr, relerr = number_error(a, b)
                if abserr > self.atol + self.rtol * max(abs(a), abs(b)):
                    errors.append(f"line {lineno}: {a!r} vs {b!r} (abs err {abserr:.3g}, rel err {relerr:.3g})"
                                  f"\n  expected: {eline}\n  actual:   {aline}")
        return errors


def measure(expected, actual):
    """
    Return a summary of the differences between two outputs: whether they are
    identical, and the largest relative and absolute differences of their numbers
    (None if the text or the number of lines or numbers differs).
    """
    if expected == actual:
        return "identical"
    elines = expected.splitlines()
    alines = actual.splitlines()
    if len(elines) != len(alines):
        return "line count differs"
    maxrel = maxabs = 0.0
    ndiff = 0
    for eline, aline in zip(elines, alines):
        eskel, enums = split_numbers(eline)
        askel, anums = split_numbers(aline)
        if eskel != askel or len(enums) != len(anums):
            return "text differs"
        for a, b in zip(enums, anums):
            abserr, relerr = number_error(a, b)
            if abserr > 0.0:
                ndiff += 1
                maxabs = max(maxabs, abserr)
                maxrel = max(maxrel, relerr)
    return f"{ndiff} numbers differ, max rel {maxrel:.2e}, max abs {maxabs:.2e}"


# ---------------------------------------------------------------------------
# test definitions
# ---------------------------------------------------------------------------

EXACT = Exact()


@dataclass
class Test:
    """
    One run of an example program.

    name     test name (the stem of the stdout reference file)
    prog     example program in the examples directory
    args     command line arguments; "{exdir}" is replaced by the reference directory
    stdout   reference file for standard output
    stderr   reference file for standard error (None: standard error must be empty)
    stdin    file in the reference directory read on standard input
    env      environment variables set for the run
    setup    example programs run first, in the same working directory
    compare  comparison of standard output (EXACT, or Numeric(rtol, atol))
    compare_err  comparison of standard error (None: same as compare)
    skip     reason for not running the test (None: run it)
    """

    name: str
    prog: str
    stdout: str
    args: tuple = ()
    stderr: str = None
    stdin: str = None
    env: dict = field(default_factory=dict)
    setup: tuple = ()
    compare: object = EXACT
    compare_err: object = None
    skip: str = None

    def comparison(self, label):
        """comparison used for the given output stream"""
        if label == "stderr" and self.compare_err is not None:
            return self.compare_err
        return self.compare


def T(prog, stdout, stderr=None, args=(), **kw):
    """Shorthand for a test named after its stdout reference file."""
    return Test(name=Path(stdout).stem, prog=prog, stdout=stdout, stderr=stderr, args=tuple(args), **kw)


# Comparison policy
#
# Outputs are compared exactly, unless their numbers depend on the floating-point
# environment. To find these, the examples were run with GSL built with
# "-O2 -march=native", which lets the compiler contract a*b+c into fused
# multiply-adds (FMA), as many compilers and platforms do by default. Building
# with -O0 or -O2 gives identical outputs.
#
# - Outputs which did not change under FMA are compared exactly.
#
# - Outputs in which FMA changes only rounding-level values (error estimates,
#   values near zero, last digits) or results of ill-conditioned problems are
#   compared numerically: rtol is the smallest power of 10 which is at least 10
#   times the largest relative difference observed, and atol likewise for values
#   near zero.
#
# - Outputs in which FMA changes a discrete choice of an adaptive algorithm (the
#   selected step size, the number of function evaluations) are compared exactly
#   and marked "platform-sensitive": a tolerance large enough to accept the
#   change would hide real changes. Their reference files may need to be
#   regenerated on platforms which use FMA.

TESTS = [
    T("blas", "blas.txt"),
    T("bspline_deriv", "bspline_deriv.txt"),
    T("bspline_gram", "bspline_gram1.txt"),
    T("bspline_gram", "bspline_gram2.txt", args=["0", "0.6"]),
    T("bspline_interp", "bspline_interp.txt"),
    T("bspline_lsbreak", "bspline_lsbreak.txt", "bspline_lsbreak.err"),
    T("bspline_lsend", "bspline_lsend.txt", "bspline_lsend.err"),
    T("bspline_lsorder", "bspline_lsorder.txt"),
    T("bspline_per", "bspline_per.txt", "bspline_per.err"),
    T("bspline_proj", "bspline_proj.txt"),
    T("cblas", "cblas.txt"),
    T("cdf", "cdf.txt"),
    T("cheb", "cheb.txt"),
    T("combination", "combination.txt"),
    T("const", "const.txt"),
    # error estimates at the rounding level (up to 6e-15 absolute)
    T("diff", "diff.txt", compare=Numeric(rtol=1e-12, atol=1e-13)),
    T("diff_d", "diff_d.txt", compare=Numeric(rtol=1e-12, atol=1e-14)),
    # platform-sensitive: under FMA the adaptive step selection chooses other steps
    T("diff_g", "diff_g.txt"),
    T("diff_j", "diff_j.txt"),
    T("diff_exp", "diff_exp.txt"),
    # Romberg error estimates: up to 2.5e-11 relative, 1.8e-14 absolute
    T("diff_diag", "diff_diag_x1.txt", "diff_diag_x1_output.txt", args=["1"],
      compare=Numeric(rtol=1e-9, atol=1e-12)),
    T("diff_diag", "diff_diag_x0.txt", "diff_diag_x0_output.txt", args=["0"]),
    T("dwt", "dwt.txt", args=["{exdir}/ecg.dat"]),
    T("eigen", "eigen.txt"),
    T("eigen_nonsymm", "eigen_nonsymm.txt"),
    # rounding noise in values near zero (up to 3e-16 / 1.1e-15 absolute)
    T("fft", "fft.txt", compare=Numeric(rtol=1e-12, atol=1e-14)),
    T("fftmr", "fftmr.txt", compare=Numeric(rtol=1e-12, atol=1e-13)),
    T("fftreal", "fftreal.txt"),
    T("filt_edge", "filt_edge.txt"),
    # ill-conditioned regularized fits: up to 4.3e-7 / 6.8e-7 relative
    T("fitreg", "fitreg.txt", "fitreg.err", compare=Numeric(rtol=1e-5)),
    T("fitreg2", "fitreg2.txt", "fitreg2.err", compare=Numeric(rtol=1e-5)),
    T("fitting", "fitting.txt"),
    T("fitting2", "fitting2.txt", args=["19"], stdin="exp.dat"),
    # last digits: up to 1.3e-13 / 9.1e-13 relative
    T("gaussfilt", "gaussfilt.txt", compare=Numeric(rtol=1e-11)),
    T("gaussfilt2", "gaussfilt2.txt", compare=Numeric(rtol=1e-11, atol=1e-15)),
    T("histogram2d", "histogram2d.txt"),
    T("ieee", "ieee.txt"),
    T("ieeeround", "ieeeround.txt"),
    T("impulse", "impulse.txt"),
    # results up to 3e-14 relative; error estimates at the rounding level (1.3e-13 absolute)
    T("integration", "integration.txt", compare=Numeric(rtol=1e-12, atol=1e-12)),
    T("integration2", "integration2a.txt", args=["10", "5"], compare=Numeric(rtol=1e-12, atol=1e-12)),
    T("integration2", "integration2b.txt", args=["10", "6"], compare=Numeric(rtol=1e-12, atol=1e-12)),
    T("interp", "interp.txt"),
    T("interp2d", "interp2d.txt"),
    T("interp_compare", "interp_compare.txt"),
    T("interpp", "interpp.txt"),
    # last digit of a Bessel function value (3e-16 relative)
    T("intro", "intro.txt", compare=Numeric(rtol=1e-14)),
    # ill-conditioned least squares (condition number 2.5e11): up to 2.1e-5 relative
    T("largefit", "largefit.txt", "largefit.err", compare=Numeric(rtol=1e-3)),
    T("largefit", "largefit2.txt", "largefit2.err", args=["1e-5"], compare=Numeric(rtol=1e-3)),
    T("linalglu", "linalglu.txt"),
    T("matrix", "matrix.txt", "matrix.err", skip="disabled in the original test script"),
    T("matrixw", "matrixw.txt"),
    T("min", "min.txt"),
    T("monte", "monte.txt"),
    T("movstat1", "movstat1.txt"),
    T("movstat2", "movstat2.txt"),
    T("movstat3", "movstat3.txt"),
    T("multiset", "multiset.txt"),
    # platform-sensitive: under FMA the solver needs one more function evaluation
    T("nlfit", "nlfit.txt", "nlfit.err"),
    # final cost of order 1e-18 varies up to 6.4e-8 relative
    T("nlfit2", "nlfit2.txt", "nlfit2.err", compare_err=Numeric(rtol=1e-6)),
    T("nlfit2b", "nlfit2b.txt", "nlfit2b.err"),
    # final cond(J) up to 5e8 varies up to 2.1e-3 relative; rtol is kept at 1e-2
    # so that a change of 1 in the iteration and evaluation counts (< 100) is detected
    T("nlfit3", "nlfit3.txt", "nlfit3.err", compare_err=Numeric(rtol=1e-2)),
    T("nlfit4", "nlfit4.txt", "nlfit4.err", skip="output contains timings"),
    T("ode-initval", "ode-initval.txt"),
    T("permseq", "permseq.txt"),
    T("permshuffle", "permshuffle.txt"),
    # residual at the rounding level (4.3e-12 absolute)
    T("poisson", "poisson.txt", "poisson.err", compare_err=Numeric(rtol=1e-12, atol=1e-10)),
    # last digits of the roots (1.2e-15 relative)
    T("polyroots", "polyroots.txt", compare=Numeric(rtol=1e-14)),
    T("qrng", "qrng.txt"),
    T("randpoisson", "randpoisson.txt"),
    T("randpoisson", "randpoisson2.txt", "randpoisson2.err", env={"GSL_RNG_SEED": "123"}),
    T("randwalk", "randwalk.txt"),
    T("rng", "rng.txt"),
    T("rngunif", "rngunif.txt"),
    T("rngunif", "rngunif2.txt", "rngunif2.err", env={"GSL_RNG_TYPE": "mrg", "GSL_RNG_SEED": "123"}),
    T("robfit", "robfit.txt", args=["100"]),
    T("rootnewt", "rootnewt.txt"),
    T("roots", "roots.txt"),
    T("rquantile", "rquantile.txt"),
    T("rstat", "rstat.txt"),
    T("siman", "siman.txt"),
    T("sortsmall", "sortsmall.txt"),
    # last digit of a Bessel function value (3e-16 relative)
    T("specfun", "specfun.txt", compare=Numeric(rtol=1e-14)),
    T("specfun_e", "specfun_e.txt", compare=Numeric(rtol=1e-14)),
    T("specfun2", "specfun2.txt"),
    T("spmatrix", "spmatrix.txt"),
    T("stat", "stat.txt"),
    T("statsort", "statsort.txt"),
    # series acceleration: sum up to 6.5e-11 relative, error estimates 1.1e-10 absolute
    T("sum", "sum.txt", compare=Numeric(rtol=1e-9, atol=1e-9)),
    T("vectorview", "vectorview.txt"),
    T("ntupler", "ntuple.txt", setup=["ntuplew"]),
    T("vectorr", "vectorr.txt", setup=["vectorw"]),
]


# ---------------------------------------------------------------------------
# running tests
# ---------------------------------------------------------------------------

# environment variables which change the output of the examples
CLEAN_ENV = ("GSL_RNG_TYPE", "GSL_RNG_SEED", "GSL_IEEE_MODE")


@dataclass
class Outcome:
    status: str               # PASS, FAIL, SKIP
    messages: list = field(default_factory=list)
    stdout: str = ""
    stderr: str = ""


def run_test(test, exdir, refdir, timeout):
    """Run one test in a temporary working directory and return its outcome."""
    prog = exdir / test.prog
    if not os.access(prog, os.X_OK):
        return Outcome("FAIL", [f"{prog} not found; build the examples with 'make -C doc/examples check'"])

    env = {k: v for k, v in os.environ.items() if k not in CLEAN_ENV}
    env.update(test.env)

    with tempfile.TemporaryDirectory(prefix="gsl_example_") as cwd:
        for name in test.setup:
            setup = subprocess.run([str(exdir / name)], cwd=cwd, env=env, capture_output=True,
                                   text=True, timeout=timeout)
            if setup.returncode != 0:
                return Outcome("FAIL", [f"setup program {name} failed with status {setup.returncode}"])

        args = [a.format(exdir=refdir) for a in test.args]
        stdin = open(refdir / test.stdin) if test.stdin else subprocess.DEVNULL
        try:
            proc = subprocess.run([str(prog)] + args, cwd=cwd, env=env, stdin=stdin,
                                  capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return Outcome("FAIL", [f"timed out after {timeout} s"])
        finally:
            if test.stdin:
                stdin.close()

    outcome = Outcome("PASS", stdout=proc.stdout, stderr=proc.stderr)

    if proc.returncode != 0:
        outcome.messages.append(f"exit status {proc.returncode}")

    for label, ref, actual in (("stdout", test.stdout, proc.stdout), ("stderr", test.stderr, proc.stderr)):
        if ref is None:
            if actual:
                outcome.messages.append(f"{label}: unexpected output (no reference file)")
            continue
        expected = (refdir / ref).read_text()
        comparison = test.comparison(label)
        errors = comparison.compare(expected, actual)
        if errors:
            outcome.messages.append(f"{label} differs from {ref} ({comparison.describe()}):")
            outcome.messages.extend("  " + e for e in errors[:10])
            if len(errors) > 10:
                outcome.messages.append(f"  ... {len(errors) - 10} more differences")

    if outcome.messages:
        outcome.status = "FAIL"

    return outcome


def main():
    root = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(description="Test the output of the GSL example programs.")
    parser.add_argument("-k", "--filter", action="append", default=[],
                        help="run only tests whose name matches this glob or substring (repeatable)")
    parser.add_argument("-v", "--verbose", action="store_true", help="show passing tests too")
    parser.add_argument("--exdir", type=Path, default=root / "doc" / "examples",
                        help="directory of the example programs (default: doc/examples)")
    parser.add_argument("--refdir", type=Path, default=root / "doc" / "examples",
                        help="directory of the reference and data files (default: doc/examples)")
    parser.add_argument("--timeout", type=float, default=120.0, help="time limit per test in seconds")
    parser.add_argument("--list", action="store_true", help="list the tests and their comparison modes")
    parser.add_argument("--measure", action="store_true",
                        help="report the largest numeric differences from the reference files")
    parser.add_argument("--update", action="store_true",
                        help="overwrite the reference files of the selected tests with the current output")
    opts = parser.parse_args()

    exdir = opts.exdir.resolve()
    refdir = opts.refdir.resolve()

    def selected(test):
        if not opts.filter:
            return True
        return any(fnmatch.fnmatch(test.name, f) or f in test.name for f in opts.filter)

    tests = [t for t in TESTS if selected(t)]

    if opts.list:
        for t in tests:
            if t.skip:
                mode = f"SKIP ({t.skip})"
            else:
                mode = t.compare.describe()
                if t.stderr is not None and t.compare_err is not None:
                    mode += f"; stderr {t.compare_err.describe()}"
            print(f"{t.name:24s} {mode}")
        return 0

    if opts.update and not opts.filter:
        parser.error("--update requires -k to select the tests whose reference files are overwritten")

    counts = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for t in tests:
        if t.skip:
            counts["SKIP"] += 1
            if opts.verbose:
                print(f"SKIP  {t.name}: {t.skip}")
            continue

        outcome = run_test(t, exdir, refdir, opts.timeout)

        if opts.measure:
            for label, ref, actual in (("stdout", t.stdout, outcome.stdout), ("stderr", t.stderr, outcome.stderr)):
                if ref is not None and (refdir / ref).exists():
                    print(f"{t.name:24s} {label}: {measure((refdir / ref).read_text(), actual)}")
            continue

        if opts.update:
            if outcome.stdout or outcome.stderr:
                (refdir / t.stdout).write_text(outcome.stdout)
                if t.stderr is not None:
                    (refdir / t.stderr).write_text(outcome.stderr)
                print(f"UPDATED {t.name}")
            else:
                print(f"NOT UPDATED {t.name}: no output ({'; '.join(outcome.messages)})")
            continue

        counts[outcome.status] += 1
        if outcome.status == "FAIL" or opts.verbose:
            print(f"{outcome.status}  {t.name}")
            for m in outcome.messages:
                print("      " + m)

    if opts.measure or opts.update:
        return 0

    ntot = counts["PASS"] + counts["FAIL"]
    print(f"\n{counts['PASS']}/{ntot} passed, {counts['FAIL']} failed, {counts['SKIP']} skipped")
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
