"""Pure-filesystem tests: no Meep, SSH, account or Git mutation."""
import importlib.util
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('export_k', ROOT / 'export_k.py')
export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(export)


class ExportTests(unittest.TestCase):
    def fixture(self, root):
        (root / 'ti2d').mkdir(parents=True)
        (root / 'ti2d/example.py').write_text('x = 1\n')
        output = root / 'runs/boundary_20260924_K/flat'
        output.mkdir(parents=True)
        raw = output / 'field.npz'
        raw.write_bytes(b'fixture-binary-evidence')
        (output / 'solver.log').write_text('fixture log\n')
        (output / 'result.json').write_text(json.dumps({'status': 'QUALIFIED', 'evidence_sha256': {'field.npz': export.sha(raw)}}))
        cache = root / 'cache/references/example'
        cache.mkdir(parents=True)
        (cache / 'reference.npz').write_bytes(b'reference-fixture')
        (root.parent / 'budget_ledger.json').write_text('{"version": 1, "intervals": []}\n')
        return output

    def test_complete_hash_checked_immutable_generations(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'K'
            self.fixture(root)
            with patch.object(export, 'ROOT', root):
                first = export.main(label='flat')
                second = export.main(final=True)
            self.assertNotEqual(first['generation'], second['generation'])
            self.assertTrue((root / 'exports' / first['raw_file']).exists())
            self.assertEqual(export.sha(root / 'exports' / second['raw_file']), second['raw_sha256'])
            self.assertEqual(json.loads((root / 'exports/export_ready.json').read_text()), second)
            with tarfile.open(root / 'exports' / second['raw_file']) as archive:
                names = archive.getnames()
                self.assertIn('boundary_absorber_K/cache/references/example/reference.npz', names)
                self.assertIn('boundary_absorber_K/runs/boundary_20260924_K/flat/field.npz', names)
            with tarfile.open(root / 'exports' / second['compact_file']) as archive:
                names = archive.getnames()
                self.assertFalse(any(name.endswith(('.npz', '.log', '.patch')) for name in names))
                self.assertIn('meep_screen_archive/boundary_absorber_K/K_FILES.json', names)
                self.assertNotIn('meep_screen_archive/FILES.json', names)

    def test_bad_qualified_hash_never_publishes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'K'
            output = self.fixture(root)
            (output / 'field.npz').write_bytes(b'changed')
            with patch.object(export, 'ROOT', root), self.assertRaisesRegex(ValueError, 'HASH_FAILED'):
                export.main()
            self.assertFalse((root / 'exports/export_ready.json').exists())

    def test_secret_in_text_blocks_raw_and_compact(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'K'
            self.fixture(root)
            (root / 'ti2d/example.py').write_text('github_pat_' + 'A' * 30)
            with patch.object(export, 'ROOT', root), self.assertRaisesRegex(ValueError, 'SECRET_SCAN'):
                export.main()
            self.assertFalse((root / 'exports/export_ready.json').exists())

    def test_symlink_does_not_follow_outside_root(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'K'
            self.fixture(root)
            (root / 'ti2d/escape').symlink_to(root.parent / 'budget_ledger.json')
            with patch.object(export, 'ROOT', root), self.assertRaisesRegex(ValueError, 'SYMLINK'):
                export.main()


if __name__ == '__main__':
    unittest.main()
