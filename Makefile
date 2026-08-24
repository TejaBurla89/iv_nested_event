.PHONY: pdf simulation test

pdf:
	pandoc guide.md \
		--from=markdown+tex_math_single_backslash+raw_tex \
		--standalone --toc --toc-depth=2 \
		--include-in-header=latex/pandoc-header.tex \
		-V geometry:margin=0.82in -V fontsize=11pt \
		-o latex/guide.tex
	xelatex -no-pdf -interaction=nonstopmode -halt-on-error \
		-output-directory=latex latex/guide.tex
	xelatex -no-pdf -interaction=nonstopmode -halt-on-error \
		-output-directory=latex latex/guide.tex
	xelatex -no-pdf -interaction=nonstopmode -halt-on-error \
		-output-directory=latex latex/guide.tex
	xdvipdfmx -o "$${TMPDIR:-/tmp}/iv_nested_event_guide.pdf" latex/guide.xdv
	cp "$${TMPDIR:-/tmp}/iv_nested_event_guide.pdf" guide.pdf

simulation:
	python3 simulation/figure3_monte_carlo.py

test:
	python3 -m unittest discover -s implementation -p 'test_*.py'
