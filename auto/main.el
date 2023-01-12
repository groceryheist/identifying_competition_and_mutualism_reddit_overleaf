(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "letterpaper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("url" "hyphens") ("hyphenat" "htt")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "aaai21"
    "times"
    "helvet"
    "courier"
    "url"
    "graphicx"
    "natbib"
    "caption"
    "color"
    "xcolor"
    "tikz"
    "booktabs"
    "multicol"
    "subcaption"
    "hyphenat"
    "amsthm"
    "amsmath"
    "commath"
    "mathtools")
   (TeX-add-symbols
    "BibTeX"
    "year"
    "UrlFont"
    "citepos"
    "oldnorm"
    "norm"
    "Slash")
   (LaTeX-add-labels
    "sec:intro"
    "sec:related.work"
    "sec:rdp"
    "sec:ecology_background"
    "sec:community_ecology"
    "sec:methods"
    "methods:density"
    "eq:user.frequency"
    "eq:user.frequency.svd"
    "eq:user.overlap"
    "eq:user.overlap.density"
    "sec:reg.H1"
    "eq:M1"
    "sec:clustering"
    "sec:var"
    "eq:var1"
    "sec:characterizing.ecological.communities"
    "eq:irf"
    "eq:average.interaction"
    "eq:average.absolute.interaction"
    "sec:mes.forecasting"
    "sec:results"
    "sec:res:studyA"
    "fig:density"
    "sec:res.characterizing"
    "fig:commense.x.abs.commense"
    "sec:case.studies"
    "fig:mut.network"
    "fig:comp.network"
    "fig:mixed.network"
    "fig:void.network"
    "fig:networks"
    "sec:res.forecasting"
    "sec:limitations"
    "sec:discussion")
   (LaTeX-add-bibliographies
    "ecological_models")
   (LaTeX-add-xcolor-definecolors
    "c77a1d2"
    "bf9837"
    "cc0c0c0"
    "mycomp"
    "mymut")
   (LaTeX-add-mathtools-DeclarePairedDelimiters
    '("norm" "")))
 :latex)

