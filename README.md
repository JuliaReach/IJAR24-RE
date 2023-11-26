# Verified propagation of imprecise probabilities in non-linear ODEs

This is the repeatability evaluation (RE) package for the article
[Verified propagation of imprecise probabilities in non-linear ODEs](https://doi.org/10.1016/j.ijar.2023.109044).
If you use this package in your work, please cite it using the metadata [here](CITATION.bib) or below.

<details>
<summary>Click to see the BibTeX entry.</summary>

```bibtex
@article{GrayFSFB24,
  author       = {Ander Gray and
                  Marcelo Forets and
                  Christian Schilling and
                  Scott Ferson and
                  Luis Benet},
  title        = {Verified propagation of imprecise probabilities in non-linear {ODE}s},
  journal      = {Int. J. Approx. Reason.},
  year         = {2024},
  volume       = {164},
  url          = {https://doi.org/10.1016/j.ijar.2023.109044},
  doi          = {10.1016/j.ijar.2023.109044},
}
```
</details>

## Generating the plots from the article

First install a Julia compiler following the instructions [here](http://julialang.org/downloads).
Once installed, download the contents of this repository and execute the following command:

```julia
$ julia --project examples/generate_plots.jl
```

The experiments are divided into subfolders in the `examples` folder.
The plots are generated in a `plots` subfolder, e.g., `examples/UnivariateOscillator/plots`.
It takes ca. 2 hours to run all experiments.
