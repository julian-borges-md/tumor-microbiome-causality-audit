# RO-2026-008 | Cross-sectional intratumoral microbial abundance
#
#   make verify    re-derive every asserted number from committed data
#   make poster    rebuild the 44x44 print poster and social cards
#   make clean     remove generated artefacts
#
# `verify` requires no network and no raw data download. It is the check that
# runs in CI. The WP0 simulation and real-data targets require the TCMA
# download described in docs/RUNBOOK.md.

PY := python3

.PHONY: install verify poster clean

install:
	$(PY) -m pip install -r requirements.txt

verify:
	$(PY) src/wp1_cross_consistency.py --datadir results --outdir results --figdir figures

poster:
	cd figures/poster && $(PY) poster_print.py && $(PY) social_cards.py

clean:
	rm -rf src/__pycache__ figures/poster/__pycache__
