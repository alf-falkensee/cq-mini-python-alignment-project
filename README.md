# cq-mini-python-alignment-project

This is a mini-python project that has been implemented in the context
of my "[CQ Beratung](https://www.cq-bildung.de)" further education in
bioinformatics.

The project consists in the extension of an input dna or protein
alignment to include a quality sequence. For instance, the alignment

```
AACTG_GTCAT
AGTCAA_CTGA
```

will be extended by a quality line such that it will be displayed to look like the following

```
AACTG_GTCAT
|::::  :::.
AGTCAA_CTGA
```

More information is given after cloning of the project using:

```
> python show_alignment.py --help
> python -i
>>> import show_alignment
>>> help(show_alignment)
```
