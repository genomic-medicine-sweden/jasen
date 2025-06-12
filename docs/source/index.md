![Jasen logo](_static/logo.png)

*Just Another System for Epityping using NGSs data*

Jasen is a pipeline for performing antibiotic resistance prediction and epityping of whole genome sequenced infectious bacteria.

# Get started

- Get an overivew of the pipeline and the different bacterial species it supports.
- Get started with installing Jasen.
- Analyze your first sample.

# Report and visualisation

* [Bonsai](https://github.com/Clinical-Genomics-Lund/cgviz): Visualises jasen outputs.
* [graptetree](https://github.com/achtman-lab/GrapeTree): Visualise phylogenetic relationship using cgmlst data.

# Tips

* Always run the latest versions of the bioinformatical software.
* Verify you have execution permission for jasens `*.sif` images.
* Old Singularity versions may sporadically produce the error `FATAL: could not open image jasen/container/*.sif: image format not recognized!`

# License

The Jasen documentation is distributed under
[GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html).


```{toctree}
:hidden:
:caption: Get started
:maxdepth: 1

overview
install
update
usage
issues
```

```{toctree}
:hidden:
:caption: Workflows
:maxdepth: 1

workflows/staphylococcus_aureus
workflows/escherichia_coli
workflows/mycobacterium_tuberculosis
workflows/klebsiella_pneumoniae
workflows/streptococcus_pyogenes
workflows/streptococcus
```

```{toctree}
:hidden:
:caption: Development
:maxdepth: 1

contributing
changelog
```