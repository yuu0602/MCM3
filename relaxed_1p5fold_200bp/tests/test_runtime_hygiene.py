"""Non-destructive syntax and R-device smoke tests; no pipeline data are used."""

import ast
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get("TEST_RSCRIPT") or shutil.which("Rscript")


def r_sources():
    for path in (ROOT / "pipeline").rglob("*.R"):
        if not path.name.startswith("._"):
            yield str(path.relative_to(ROOT)), path.read_text()
    for path in (ROOT / "pipeline").rglob("*.py"):
        if path.name.startswith("._"):
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "r_script"
                for target in node.targets
            ):
                yield str(path.relative_to(ROOT)), ast.literal_eval(node.value)


class RuntimeHygiene(unittest.TestCase):
    def test_python_syntax_and_r_paths(self):
        for path in ROOT.rglob("*"):
            if path.name.startswith("._") or path.suffix not in (".py", ".sh"):
                continue
            source = path.read_text()
            self.assertNotIn("/usr/local/bin/" + "Rscript", source, str(path))
            if path.suffix == ".py":
                ast.parse(source, filename=str(path))

    @unittest.skipUnless(RSCRIPT, "Rscript not available")
    def test_r_syntax_and_device_cleanup(self):
        sources = list(r_sources())
        self.assertEqual(len(sources), 5)
        for name, source in sources:
            with self.subTest(script=name):
                proc = subprocess.run(
                    [RSCRIPT, "--vanilla", "-e", 'parse(text=readLines(file("stdin")))'],
                    input=source, text=True, capture_output=True,
                )
                self.assertEqual(proc.returncode, 0, proc.stderr)
                marker = "  grDevices::pdf(file = NULL)\n"
                guard = source[source.index("local({"):source.index(marker) + len(marker)]
                for fail in (False, True):
                    with tempfile.TemporaryDirectory(prefix="mcm3-r-devices-") as directory:
                        # Only disposable fixtures are created or removed by this test.
                        body = '''
plot.new()
png("intended.png", width=600, height=600)
plot(1:3)
dev.off()
writeLines("default artifact fixture", "Rplots.pdf")
writeLines("resource fork fixture", "._Rplots.pdf")
'''
                        if fail:
                            body += 'stop("intentional smoke-test failure")\n'
                        code = 'before <- getOption("device")\ntryCatch(\n' + guard + body
                        code += '}), error=function(e) NULL)\n'
                        code += 'stopifnot(is.null(dev.list()), identical(before, getOption("device")))\n'
                        proc = subprocess.run(
                            [RSCRIPT, "--vanilla", "-"], input=code,
                            text=True, capture_output=True, cwd=directory,
                        )
                        self.assertEqual(proc.returncode, 0, proc.stderr)
                        self.assertGreater((Path(directory) / "intended.png").stat().st_size, 0)
                        self.assertFalse((Path(directory) / "Rplots.pdf").exists())
                        self.assertFalse((Path(directory) / "._Rplots.pdf").exists())

    def test_stage06_resolves_r_against_child_path(self):
        path = ROOT / "pipeline/regulatory_target/06_define_regulatory_targets.py"
        spec = importlib.util.spec_from_file_location("stage06_runtime_test", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory(prefix="mcm3-path-") as directory:
            with patch.object(module, "REGULATORY_ROOT", Path(directory)), \
                 patch.object(module.sys, "executable", "/mock-conda/bin/python"), \
                 patch.object(module.shutil, "which", return_value="/mock-conda/bin/Rscript") as which, \
                 patch.object(module.subprocess, "run") as run:
                module.run_r("render_direction.R")
                self.assertTrue(which.call_args.kwargs["path"].startswith("/mock-conda/bin" + os.pathsep))
                self.assertEqual(run.call_args.args[0][0], "/mock-conda/bin/Rscript")
                self.assertEqual(which.call_args.kwargs["path"], run.call_args.kwargs["env"]["PATH"])


if __name__ == "__main__":
    unittest.main()
