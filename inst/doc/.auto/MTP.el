(TeX-add-style-hook "MTP"
 (lambda ()
    (LaTeX-add-bibliographies
     "xing")
    (LaTeX-add-environments
     "theorem"
     "procedure")
    (LaTeX-add-labels
     "anal:mult:multtest"
     "anal:mult:s:intro"
     "anal:mult:s:methods"
     "anal:mult:s:framework"
     "anal:mult:s:nullDistn"
     "anal:mult:s:software"
     "anal:mult:s:disc")
    (TeX-add-symbols
     '("Rclass" 1)
     '("Robject" 1)
     '("Rpackage" 1)
     "RR"
     "ZZ"
     "NN")
    (TeX-run-style-hooks
     "natbib"
     "authoryear"
     "round"
     "comment"
     "color"
     "amsmath"
     "hyperref"
     "amsfonts"
     "Sweave"
     "graphicx"
     "latex2e"
     "art11"
     "article"
     "11pt")))

