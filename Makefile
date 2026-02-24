.PHONY: docs-setup docs docs-serve docs-clean

docs-setup:
	julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

docs: docs-setup
	julia --project=docs docs/make.jl

docs-serve:
	cd docs/build && python3 -m http.server 8000

docs-clean:
	rm -rf docs/build
