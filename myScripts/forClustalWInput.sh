#!/bin/bash
find -type f -name "*fsta" -exec cat {} \; | dos2unix | tr -d [\\n] | sed 's/>/\'$'\n>/g' | sed 's/ab1/&\n/g' > inputToClustalW.fsta

