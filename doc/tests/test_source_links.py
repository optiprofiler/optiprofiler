"""Exercise the actual Sphinx config using tiny isolated source checkouts.

No OptiProfiler numerical dependency or Sphinx installation is needed. Git
commits below belong only to disposable fixtures, never to the project tree.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


CONFIG = Path(__file__).resolve().parents[1] / 'source' / 'conf.py'
CORE_ENV = 'OPTIPROFILER_DOCS_SOURCE_REF'
S2MPJ_ENV = 'OPTIPROFILER_DOCS_S2MPJ_REF'


class SourceLinksTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='op-doc-links-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / 'checkout'
        self.package = self.root / 'python' / 'optiprofiler'
        self.provider = self.package / 'problem_libs' / 's2mpj'
        self.provider.mkdir(parents=True)
        (self.root / 'doc' / 'source').mkdir(parents=True)
        shutil.copy2(CONFIG, self.root / 'doc' / 'source' / 'conf.py')
        (self.package / '__init__.py').write_text("__version__ = '1.3.4'\n")
        (self.package / 'core.py').write_text('def example():\n    return 1\n')
        (self.provider / 'example.py').write_text('def example():\n    return 2\n')
        (self.root / '.gitignore').write_text('python/optiprofiler/problem_libs/s2mpj/\n')
        self.core_sha = self.init_repo(self.root)
        self.provider_sha = self.init_repo(self.provider)

    @staticmethod
    def git(root, *args):
        env = os.environ.copy()
        for key in ('GIT_DIR', 'GIT_WORK_TREE', 'GIT_COMMON_DIR', 'GIT_INDEX_FILE'):
            env.pop(key, None)
        return subprocess.check_output(
            ['git', '-C', str(root), *args], env=env, text=True,
            stderr=subprocess.DEVNULL).strip()

    def init_repo(self, root):
        self.git(root, 'init', '-q')
        self.git(root, 'add', '.')
        self.git(root, '-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
                 '-c', 'commit.gpgsign=false', 'commit', '-qm', 'fixture')
        return self.git(root, 'rev-parse', 'HEAD')

    def load_config(self, env_delta=None, root=None, preload=False):
        root = self.root if root is None else root
        env = os.environ.copy()
        for key in list(env):
            if key.startswith('OPTIPROFILER_DOCS_') or key == 'PYTHONPATH':
                env.pop(key)
        env.update(env_delta or {})
        script = f'''import importlib, json, runpy, sys
{'import optiprofiler' if preload else ''}
conf = runpy.run_path({str(root / 'doc' / 'source' / 'conf.py')!r})
import optiprofiler
importlib.import_module('optiprofiler.core')
importlib.import_module('optiprofiler.problem_libs.s2mpj.example')
resolve = conf['linkcode_resolve']
print(json.dumps(dict(
    package=str(optiprofiler.__file__), release=conf['release'],
    core=resolve('py', {{'module':'optiprofiler.core', 'fullname':'example'}}),
    provider=resolve('py', {{'module':'optiprofiler.problem_libs.s2mpj.example', 'fullname':'example'}}),
    external=resolve('py', {{'module':'optiprofiler.core', 'fullname':'external'}}))))
'''
        # The docs config must not depend on the shell's current directory.
        return subprocess.run([sys.executable, '-c', script], cwd=self.base,
                              env=env, text=True, capture_output=True)

    def output(self, **kwargs):
        result = self.load_config(**kwargs)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout.strip().splitlines()[-1])

    def test_core_and_provider_link_their_actual_revisions(self):
        output = self.output(env_delta={'PYTHONPATH': str(self.root / 'python')})
        self.assertEqual(output['core'], 'https://github.com/optiprofiler/optiprofiler/'
                         f'blob/{self.core_sha}/python/optiprofiler/core.py#L1-L2')
        self.assertEqual(output['provider'], 'https://github.com/optiprofiler/s2mpj_python/'
                         f'blob/{self.provider_sha}/example.py#L1-L2')
        self.assertEqual(output['release'], '1.3.4')

    def test_source_path_precedes_an_installed_package(self):
        installed = self.base / 'installed'
        package = installed / 'optiprofiler'
        package.mkdir(parents=True)
        (package / '__init__.py').write_text("__version__='9.9.9'\n")
        output = self.output(env_delta={'PYTHONPATH': str(installed)})
        self.assertEqual(Path(output['package']).resolve(), (self.package / '__init__.py').resolve())
        self.assertEqual(output['release'], '1.3.4')

    def test_already_imported_foreign_package_is_rejected(self):
        installed = self.base / 'installed'
        package = installed / 'optiprofiler'
        package.mkdir(parents=True)
        (package / '__init__.py').write_text("__version__='9.9.9'\n")
        result = self.load_config(env_delta={'PYTHONPATH': str(installed)}, preload=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not belong to this documentation checkout', result.stderr)

    def test_gitless_archive_has_no_guessed_links(self):
        shutil.rmtree(self.root / '.git')
        shutil.rmtree(self.provider / '.git')
        output = self.output()
        self.assertIsNone(output['core'])
        self.assertIsNone(output['provider'])

    def test_archive_does_not_borrow_an_unrelated_parent_revision(self):
        shutil.rmtree(self.root / '.git')
        shutil.rmtree(self.provider / '.git')
        self.init_repo(self.base)
        output = self.output()
        self.assertIsNone(output['core'])
        self.assertIsNone(output['provider'])

    def test_git_unavailable_has_no_guessed_links(self):
        output = self.output(env_delta={'PATH': ''})
        self.assertIsNone(output['core'])
        self.assertIsNone(output['provider'])

    def test_provider_without_git_does_not_borrow_core_revision(self):
        shutil.rmtree(self.provider / '.git')
        output = self.output()
        self.assertIn(f'/blob/{self.core_sha}/', output['core'])
        self.assertIsNone(output['provider'])

    def test_reexport_outside_checkout_has_no_source_link(self):
        with (self.package / 'core.py').open('a') as stream:
            stream.write('from pathlib import Path as external\n')
        self.assertIsNone(self.output()['external'])

    def test_gitless_explicit_full_sha_overrides(self):
        shutil.rmtree(self.root / '.git')
        shutil.rmtree(self.provider / '.git')
        output = self.output(env_delta={CORE_ENV: self.core_sha, S2MPJ_ENV: self.provider_sha})
        self.assertIn(f'/blob/{self.core_sha}/', output['core'])
        self.assertIn(f'/blob/{self.provider_sha}/', output['provider'])

    def test_override_cannot_disagree_with_checkout(self):
        for key in (CORE_ENV, S2MPJ_ENV):
            with self.subTest(key=key):
                result = self.load_config(env_delta={key: '0' * 40})
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('does not match the checked-out revision', result.stderr)

    def test_matching_overrides_are_accepted(self):
        output = self.output(env_delta={CORE_ENV: self.core_sha.upper(),
                                        S2MPJ_ENV: self.provider_sha.upper()})
        self.assertIn(f'/blob/{self.core_sha}/', output['core'])
        self.assertIn(f'/blob/{self.provider_sha}/', output['provider'])

    def test_invalid_refs_are_rejected_not_interpreted(self):
        for value in ('main', 'paper/v1.3.4-errata', 'abc1234', '../main', 'a' * 40 + '?x=1'):
            with self.subTest(value=value):
                result = self.load_config(env_delta={CORE_ENV: value})
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('full 40-character hexadecimal commit SHA', result.stderr)

    def test_git_worktree_file_is_supported(self):
        worktree = self.base / 'worktree'
        self.git(self.root, 'worktree', 'add', '--detach', '-q', str(worktree), self.core_sha)
        # Populate the ignored provider independently, as a normal submodule would be.
        target = worktree / 'python' / 'optiprofiler' / 'problem_libs' / 's2mpj'
        shutil.copytree(self.provider, target)
        output = self.output(root=worktree)
        self.assertIn(f'/blob/{self.core_sha}/', output['core'])

    def test_inherited_git_environment_cannot_change_repository_identity(self):
        foreign = self.base / 'foreign'
        foreign.mkdir()
        (foreign / 'unrelated.txt').write_text('unrelated\n')
        self.init_repo(foreign)
        output = self.output(env_delta={'GIT_DIR': str(foreign / '.git'),
                                       'GIT_WORK_TREE': str(foreign)})
        self.assertIn(f'/blob/{self.core_sha}/', output['core'])
        self.assertIn(f'/blob/{self.provider_sha}/', output['provider'])


if __name__ == '__main__':
    unittest.main()
