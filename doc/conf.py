# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sphinx_rtd_theme

project = 'GSL'
copyright = '1996-2026 The GSL Team'
author = 'The GSL Team'
release = '2.8'
title = 'GNU Scientific Library'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

primary_domain = 'c'
numfig = True

extensions = ['sphinx.ext.imgmath']

# The master toctree document.
master_doc = 'index'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'include.rst', 'specfunc-*.rst']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for LaTeX output ---------------------------------------------

my_latex_preamble = '\\DeclareMathOperator\\arccosh{arccosh} \
                     \\DeclareMathOperator\\arcsinh{arcsinh} \
                     \\DeclareMathOperator\\arctanh{arctanh} \
                     \\DeclareMathOperator\\arcsec{arcsec} \
                     \\DeclareMathOperator\\arccsc{arccsc} \
                     \\DeclareMathOperator\\arccot{arccot} \
                     \\DeclareMathOperator\\csch{csch} \
                     \\DeclareMathOperator\\sech{sech} \
                     \\DeclareMathOperator\\arcsech{arcsech} \
                     \\DeclareMathOperator\\arccsch{arccsch} \
                     \\DeclareMathOperator\\arccoth{arccoth} \
                     \\DeclareMathOperator\\erf{erf} \
                     \\DeclareMathOperator\\erfc{erfc} \
                     \\DeclareMathOperator\\sgn{sgn} \
                     \\DeclareMathOperator\\sinc{sinc} \
                     \\DeclareMathOperator\\Var{Var} \
                     \\DeclareMathOperator\\diag{diag} \
                     \\newcommand\\undermat[2]{\\makebox[0pt][l]{$\\smash{\\underbrace{\\phantom{\\begin{matrix}#2\\end{matrix}}}_{\\text{$#1$}}}$}#2}'

my_latex_authors = 'Mark Galassi \\\\ \
                    Jim Davies \\\\ \
                    James Theiler \\\\ \
                    Brian Gough \\\\ \
                    Gerard Jungman \\\\ \
                    Patrick Alken \\\\ \
                    Michael Booth \\\\ \
                    Fabrice Rossi \\\\ \
                    Rhys Ulerich'

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    'preamble': my_latex_preamble,

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'gsl-ref.tex', title,
     my_latex_authors, 'manual'),
]

imgmath_latex_preamble = my_latex_preamble

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'gsl', title,
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'gsl-ref', title,
     author, 'GSL', 'One line description of project.',
     'Miscellaneous'),
]
