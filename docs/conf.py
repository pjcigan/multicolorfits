# -- Path setup --------------------------------------------------------------
import os
import pathlib
import re
import sys
import warnings

sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# nbsphinx still passes string paths through Sphinx's doc2path; Sphinx 8 emits
# RemovedInSphinx90Warning per notebook.  Filter until upstream is Path-clean.
try:
    from sphinx.deprecation import RemovedInSphinx90Warning
except ImportError:  # pragma: no cover
    RemovedInSphinx90Warning = None  # type: ignore[misc, assignment]
if RemovedInSphinx90Warning is not None:
    warnings.filterwarnings(
        'ignore', category=RemovedInSphinx90Warning, module=r'nbsphinx(\.|$)')

# Register denim Pygments styles (same machinery as skyplothelper).
from pygments.styles import STYLES, _STYLE_NAME_TO_MODULE_MAP  # noqa: E402

STYLES['DenimDarkStyle'] = ('_pygments_denim', 'denimdark', ())
STYLES['DenimLightStyle'] = ('_pygments_denim', 'denimlight', ())
_STYLE_NAME_TO_MODULE_MAP['denimdark'] = ('_pygments_denim', 'DenimDarkStyle')
_STYLE_NAME_TO_MODULE_MAP['denimlight'] = ('_pygments_denim', 'DenimLightStyle')

# -- Project information -----------------------------------------------------
project = 'multicolorfits'
copyright = 'Phil Cigan'
author = 'Phil Cigan'

_init = (pathlib.Path(__file__).parent.parent / 'multicolorfits' / '__init__.py').read_text()
_m = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", _init)
release = _m.group(1) if _m else '3.0.0'
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    'myst_parser',
    'nbsphinx',
    'sphinx_design',
    'sphinx_copybutton',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
]

autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_rtype = False

templates_path = ['_templates']
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    '.docs_test_env',
    '**.ipynb_checkpoints',
    'tutorials/*.py',
    'README.md',
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

myst_enable_extensions = ['colon_fence', 'deflist']
myst_heading_anchors = 3

nbsphinx_execute = 'never'
nbsphinx_allow_errors = True
nbsphinx_prolog = """
{% set docname = env.doc2path(env.docname, base=None) %}
.. note::

   This page is generated from a Jupyter notebook.
   Download the ``.ipynb`` from the sources, or open the paired jupytext
   ``.py`` under ``docs/tutorials/``.
"""

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
}

# Optional deps — mock so autodoc imports succeed on RTD without extras
autodoc_mock_imports = [
    'fastapi', 'uvicorn', 'PySide6', 'reproject', 'skyplothelper', 'pysymlog',
]

# -- HTML output -------------------------------------------------------------
html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']
html_js_files = ['plot-theme.js']
html_title = 'MultiColorFits'
# Serve the agent-facing map + recipe corpus at the site root, per the llms.txt
# convention (agents fetch https://<site>/llms.txt). Both are generated from
# multicolorfits/_overview.py by scripts/make_llms_txt.py — regenerate them,
# never hand-edit.
html_extra_path = ['../llms.txt', '../llms-full.txt']
# Branding marks live in ``_static/logo/`` (regenerate with ``python docs/make_logo.py``).
# Favicon = bright RYB combo_swatch; navbar mark is composed in
# ``_templates/navbar-logo.html`` (color filled circles by default; mono outline
# variants also generated for a future swap).  html_logo intentionally omitted.
html_favicon = '_static/logo/logo_favicon.png'

html_theme_options = {
    'github_url': 'https://github.com/pjcigan/multicolorfits',
    'icon_links': [
        {
            'name': 'PyPI',
            'url': 'https://pypi.org/project/multicolorfits/',
            'icon': 'fa-brands fa-python',
        },
    ],
    'use_edit_page_button': False,
    'navbar_align': 'left',
    'pygments_light_style': 'denimlight',
    'pygments_dark_style': 'denimdark',
    'header_links_before_dropdown': 8,
    'show_toc_level': 2,
    'navigation_with_keys': True,
}

# Version switcher only on Read the Docs (needs hosted switcher.json + >1 build).
if os.environ.get('READTHEDOCS') == 'True':
    switcher_version = os.environ.get('READTHEDOCS_VERSION') or 'latest'
    html_theme_options['switcher'] = {
        'json_url': 'https://multicolorfits.readthedocs.io/en/latest/_static/switcher.json',
        'version_match': switcher_version,
    }
    html_theme_options['show_version_warning_banner'] = True
    html_theme_options['navbar_end'] = [
        'version-switcher',
        'theme-switcher',
        'navbar-icon-links',
    ]

html_context = {
    'github_user': 'pjcigan',
    'github_repo': 'multicolorfits',
    'github_version': 'master',
    'doc_path': 'docs',
}

html_sidebars = {
    '**': ['sidebar-nav-bs'],
}