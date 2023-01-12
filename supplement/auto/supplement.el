(TeX-add-style-hook
 "supplement"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("acmart" "manuscript")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "acmart"
    "acmart10"
    "color"
    "rotating"
    "tikz"
    "xcolor"
    "subfig")
   (TeX-add-symbols
    "BibTeX")
   (LaTeX-add-labels
    "mut.coefs"
    "comp.coefs"
    "mixed.coefs"
    "diet.coefs"
    "mut.irf"
    "comp.irf"
    "mixed.irf"
    "void.irf")
   (LaTeX-add-xcolor-definecolors
    "mutualism"
    "competition"))
 :latex)

