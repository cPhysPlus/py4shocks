# Import libraries
import os
import sys

# If your extensions live in 'your_module', put its path here.
sys.path.insert(0, os.path.abspath('../..'))

project = 'Your Project Name'
copyright = '2025, Yachay Tech'
author = 'B. Pinargore, W. Banda-Barragan, S. Navarrete, F. Teutloff'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = ['_build']

html_theme = 'sphinx_rtd_theme' # Good for RTD
