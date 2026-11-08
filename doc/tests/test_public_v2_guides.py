"""Keep shared v2 guides reachable in the published documentation tree."""
from pathlib import Path
import re
import unittest


SOURCE = Path(__file__).resolve().parents[1] / 'source'


class PublicV2Guides(unittest.TestCase):
    def test_shared_guides_are_in_the_user_navigation(self):
        navigation = (SOURCE / 'user' / 'index.rst').read_text(encoding='utf-8')
        for name in ('reference_targets', 'plot_options'):
            with self.subTest(guide=name):
                self.assertTrue((SOURCE / 'user' / (name + '.rst')).is_file())
                self.assertIn('    ' + name + '\n', navigation)

    def test_reference_links_do_not_target_the_excluded_dev_tree(self):
        # dev/ is deliberately excluded by conf.py. A source file existing
        # there is not enough to make a public :doc: link work on the website.
        pages = ('user/usage.rst', 'user/usage_python.rst',
                 'user/problem_libraries/custom_python.rst',
                 'user/problem_libraries/custom_matlab.rst')
        for page in pages:
            path = SOURCE / page
            links = re.findall(r':doc:`([^`]+)`', path.read_text(encoding='utf-8'))
            for link in links:
                with self.subTest(page=page, link=link):
                    target = (path.parent / (link + '.rst')).resolve()
                    self.assertTrue(target.is_file(), str(target))
                    self.assertNotIn(SOURCE / 'dev', target.parents)


if __name__ == '__main__':
    unittest.main()
