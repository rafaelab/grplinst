---
layout: default
title: References
---

# References

## Main Project Paper

- R. Alves Batista, A. Saveliev, E. M. de Gouveia Dal Pino, "Plasma effects on fast pair beams and application to the intergalactic medium", MNRAS 489 (2019) 3836. [MNRAS](https://academic.oup.com/mnras/article/489/3/3836/5558261), [arXiv:1904.13345](https://arxiv.org/abs/1904.13345)

## Implemented Model References

- P. Broderick, P. Chang, C. Pfrommer, Astrophys. J. 752 (2012) 22.
- R. Schlickeiser, D. Ibscher, M. Supsar, Astrophys. J. 758 (2012) 102.
- L. Sironi, D. Giannios, Astrophys. J. 787 (2014) 49. [arXiv:1312.4538](https://arxiv.org/abs/1312.4538)
- S. Vafin, M. Pohl, J. Niemiec, A. Bret, Astrophys. J. 865 (2018) 23. [arXiv:1807.04203](https://arxiv.org/abs/1807.04203)
- A. Bret, L. Gremillet, M. Dieckmann, Phys. Plasmas 17 (2010) 120501. [arXiv:1010.5763](https://arxiv.org/abs/1010.5763)
- M. Shalaby et al., Phys. Rev. Lett. 124 (2020) 105101. [arXiv:1907.13350](https://arxiv.org/abs/1907.13350)
- F. Miniati, A. Elyiv, Astrophys. J. 770 (2013) 54. [arXiv:1208.1761](https://arxiv.org/abs/1208.1761)

## Project Links

- CRPropa: [https://github.com/CRPropa/CRPropa3](https://github.com/CRPropa/CRPropa3)
- `grplinst`: [https://github.com/rafaelab/grplinst](https://github.com/rafaelab/grplinst)

## Repository Files Worth Reading

- Public umbrella header: [include/grplinst.h](https://github.com/rafaelab/grplinst/blob/v2/include/grplinst.h)
- Plasma-instability interface: [include/grplinst/PlasmaInstability.h](https://github.com/rafaelab/grplinst/blob/v2/include/grplinst/PlasmaInstability.h)
- Main implementation file: [src/PlasmaInstability.cc](https://github.com/rafaelab/grplinst/blob/v2/src/PlasmaInstability.cc)
- Example simulation: [examples/testPlugin.py](https://github.com/rafaelab/grplinst/blob/v2/examples/testPlugin.py)
- Python binding tests: [test/testPlugin.py](https://github.com/rafaelab/grplinst/blob/v2/test/testPlugin.py)

## Reading Order

If you are new to the codebase, a productive order is:

1. read the main project paper,
2. inspect [Models](models.html) to identify the implemented prescription you need,
3. inspect [Usage](usage.html) for the CRPropa integration pattern,
4. inspect the linked source files if you need the exact formulae.