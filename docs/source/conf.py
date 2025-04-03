import os
import sys

# If your extensions live in 'your_module', put its path here.
sys.path.insert(0, os.path.abspath('../..'))

project = 'Your Project Name'
copyright = '2024, Your Name'
author = 'Your Name'

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

#Try to import the sphinx_rtd_theme
try:
    import sphinx_rtd_theme
    print("sphinx_rtd_theme imported successfully")
except ImportError:
    print("sphinx_rtd_theme import failed")
