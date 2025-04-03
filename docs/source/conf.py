# Import librariea
import os
import sys

# Path for extensions
sys.path.insert(0, os.path.abspath('../..'))

project = 'py4shocks'
copyright = '2025, Yachay Tech University'
author = 'B. Pinargote, W. E. Banda-Barragan, S. Navarrete, F. Teutloff'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = ['_build']

html_theme = 'sphinx_rtd_theme' # Good for RTD

# Add print statements to check the python path
print("Python Path:")
for path in sys.path:
    print(path)
