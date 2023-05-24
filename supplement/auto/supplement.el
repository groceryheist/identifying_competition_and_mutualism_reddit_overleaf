(TeX-add-style-hook
 "supplement"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("acmart" "manuscript")))
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

