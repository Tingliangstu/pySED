# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import re
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pySED'
copyright = '2025, Ting Liang; Wenwu Jiang'
author = 'Ting Liang; Wenwu Jiang'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['recommonmark', 'sphinx_markdown_tables', 'nbsphinx', 'sphinx.ext.mathjax']

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
}


def _markdown_links_to_rst(text):
    def convert_link(match):
        label = match.group(1)
        url = match.group(2)
        label = re.sub(r'<sub>(.*?)</sub>', r'\1', label)
        label = label.replace('&theta;', 'theta').replace('\\&', '&')
        return f'`{label} <{url}>`_'

    return re.sub(r'\[([^\]]+)\]\(([^)]+)\)', convert_link, text)


def _sync_publications_page():
    """Create the docs publications page from the repository publication list."""
    docs_source = Path(__file__).resolve().parent
    repo_root = docs_source.parents[1]
    source = repo_root / 'publications' / 'readme.md'
    target = docs_source / 'publications.rst'
    legacy_markdown_target = docs_source / 'publications.md'

    if not source.exists():
        return

    if legacy_markdown_target.exists():
        legacy_markdown_target.unlink()

    data = source.read_bytes()
    for encoding in ('utf-8', 'gbk', 'cp1252'):
        try:
            text = data.decode(encoding)
            break
        except UnicodeDecodeError:
            text = None
    if text is None:
        text = data.decode('utf-8', errors='replace')
    lines = text.replace('\r\n', '\n').replace('\r', '\n').splitlines()
    output = []
    for line in lines:
        stripped = line.strip()
        if stripped == '# Publications using pySED':
            title = 'Publications'
            output.extend([
                title,
                '=' * len(title),
                '',
                'This page is generated from '
                '`publications/readme.md <https://github.com/Tingliangstu/pySED/tree/main/publications>`_. '
                'To update this page, edit ``publications/readme.md``, rebuild the documentation, '
                'and push the change to GitHub. ReadTheDocs will sync it during the next build.',
                '',
            ])
        elif stripped.startswith('## '):
            title = stripped[3:]
            output.extend(['', title, '-' * len(title), ''])
        elif stripped.startswith('* '):
            output.append('- ' + _markdown_links_to_rst(stripped[2:]))
        elif stripped:
            output.append(_markdown_links_to_rst(stripped))
        elif output and output[-1] != '':
            output.append('')

    target.write_text('\n'.join(output).rstrip() + '\n', encoding='utf-8')


_sync_publications_page()



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"

#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Specify the path to the static files
html_static_path = ['_static']

# Specify the logo file
html_logo = '_static/logo.png'  # Ensure 'logo.png' is placed inside the '_static' directory

# Theme options
html_theme_options = {
    'logo_only': False,           # Show project name alongside logo
    #'display_version': True,      # Display the version number
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    # Add other options as needed
}

# GitHub context for the edit link
html_context = {
    'display_github': True,  # Enable the display of the GitHub link
    'github_user': 'Tingliangstu',  # Replace with your GitHub username
    'github_repo': 'pySED',  # Replace with your repository name
    'github_version': 'main',  # Your main branch name (usually 'main' or 'master')
    'conf_py_path': '/docs/source/',  # Path to your conf.py file
}

# Custom CSS files
html_css_files = [
    'custom.css',  # Ensure 'custom.css' is placed inside the '_static' directory
]
