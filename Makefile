.PHONY: install install-fit reference website serve bench bench-ols bench-quick

install: ## Install the package
	R CMD INSTALL .

# The website covers xtx and xtx.fit (the model fits): its pages, the
# fitting benchmarks (website/benchmarks/fixed_effects, the dev harness,
# bench/compression) and the reference pages need xtx.fit installed, and
# its sources at XTX_FIT_DIR (default ../xtx.fit) for the reference pages.
XTX_FIT_DIR ?= ../xtx.fit

install-fit: ## Install xtx.fit from XTX_FIT_DIR (after xtx)
	R CMD INSTALL $(XTX_FIT_DIR)

reference: ## Convert man/*.Rd of xtx and xtx.fit to website/reference*/*.typ
	XTX_FIT_DIR=$(XTX_FIT_DIR) ./website/reference.sh

# Benchmarks never run when the site builds: the pages read saved results.
# BENCH_SIZE=full (default) saves the published results; BENCH_SIZE=quick
# runs small grids and saves them in $BENCH_OUT/website and $BENCH_OUT/ols
# (default BENCH_OUT: $TMPDIR/xtx-bench-quick), never over the
# published ones. The fixed-effects data cache is outside the repository too
# (FE_BENCH_CACHE, see website/benchmarks/fixed_effects.R).
BENCH_SIZE ?= full

bench: install install-fit ## Run website benchmarks (all, or B="crossprod qr"; BENCH_SIZE=quick|full)
	BENCH_SIZE=$(BENCH_SIZE) Rscript website/benchmarks/run.R $(B)

bench-ols: install ## Run the least squares comparison (fresh process per method; BENCH_SIZE=quick|full)
	BENCH_SIZE=$(BENCH_SIZE) ./bench/ols/run.sh

bench-quick: install install-fit ## Quick check of every published benchmark; results outside the repository
	BENCH_SIZE=quick Rscript website/benchmarks/run.R $(B)
	BENCH_SIZE=quick ./bench/ols/run.sh

# DUCKDB_R_HOME: duckdb's own default folder, named explicitly so that it
# does not print a notice about it (with its path) in the output of the pages.
# GitHub Pages serves docs/ on main: the built site is committed there
# (.nojekyll: files and folders that start with an underscore are served).
website: install install-fit reference ## Build the site into docs/ (GitHub Pages)
	DUCKDB_R_HOME=$$HOME/.duckdb calepin compile website docs
	touch docs/.nojekyll

serve: ## Serve the built site
	calepin serve docs --open
