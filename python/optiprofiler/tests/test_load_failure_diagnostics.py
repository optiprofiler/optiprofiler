"""A problem-load failure is logged with its (shortened) cause."""
import numpy as np

from optiprofiler import benchmark


def finish(fun, x0):
    x = np.ones_like(x0)
    fun(x)
    return x


def stay(fun, x0):
    fun(x0)
    return x0


def test_load_failure_log_names_the_cause(tmp_path, monkeypatch, caplog, capsys):
    monkeypatch.chdir(tmp_path)
    library = tmp_path / 'libraries' / 'failing_load'
    library.mkdir(parents=True)
    (library / 'failing_load_tools.py').write_text(
        "from optiprofiler import Problem\n"
        "def objective(x): return float(x @ x)\n"
        "def failing_load_select(options): return ['GOOD', 'BROKEN']\n"
        "def failing_load_load(name):\n"
        "    if name == 'BROKEN':\n"
        "        raise ImportError('undefined symbol _ZGVbN2v_cos (simulated)')\n"
        "    return Problem(objective, [1.0], name=name)\n",
        encoding='utf-8')
    with caplog.at_level('WARNING'):
        benchmark([finish, stay], plibs=['failing_load'], custom_problem_libs_path=str(library.parent),
                  ptype='u', mindim=1, maxdim=1, n_jobs=1, max_eval_factor=2, score_only=True,
                  draw_hist_plots='none')
    captured = capsys.readouterr()
    text = caplog.text + captured.out + captured.err
    assert 'Failed to load' in text and 'BROKEN' in text
    assert 'undefined symbol _ZGVbN2v_cos (simulated)' in text
