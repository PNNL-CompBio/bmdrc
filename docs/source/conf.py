# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'bmdrc'
copyright = '2025, David Degnan, Lisa Bramer, Lisa Truong, Robyn Tanguary, Sara Gosline, Katrina Waters'
author = 'David Degnan, Lisa Bramer, Lisa Truong, Robyn Tanguary, Sara Gosline, Katrina Waters'
release = '1.0.0'
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "sphinx.ext.coverage", "sphinx.ext.napoleon"]

autosummary_generate = True
templates_path = ['_templates']
exclude_patterns = []
autodoc_typehints = ["none"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_static_path = ['_static']