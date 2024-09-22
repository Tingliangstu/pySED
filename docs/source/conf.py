# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pySED'
copyright = '2024, Ting Liang; Wenwu Jiang'
author = 'Ting Liang; Wenwu Jiang'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['recommonmark','sphinx_markdown_tables'] 

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Specify the path to the static files
html_static_path = ['_static']

# Specify the logo file
html_logo = '_static/logo.png'  # Ensure 'logo.png' is placed inside the '_static' directory

# Theme options
html_theme_options = {
    'logo_only': False,           # Show project name alongside logo
    'display_version': True,      # Display the version number
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    # Add other options as needed
}

# Custom CSS files
html_css_files = [
    'custom.css',  # Ensure 'custom.css' is placed inside the '_static' directory
]
