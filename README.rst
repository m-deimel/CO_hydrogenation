Selectivity Trends and Role of Adsorbate-Adsorbate Interactions in CO Hydrogenation on Rhodium Catalysts
========================================================================================================

This repository contains the input files for the |kmos_paper|_ calculations together with the geometries and data used in the paper |CO_hydrogenation|_.

KMC Models
----------

The different kinetic Monte Carlo (KMC) models for |kmos_paper|_ are in the folder ``KMC_models``. To run the models clone this |kmos|_ and checkout the branch |kmos_branch|_.

    ``Rh111``

    This folder contains the python and xml files of the Rh(111) model with lateral interactions.
\

    ``Rh211``

    This folder contains the python and xml files of the Rh(211) models with and without lateral interactions.
\

    ``tools``

    This folder contains a python tool for |kmos_paper|_ to include nearest neighbors as bystanders and modify the otf_rate accordingly.

Structures
----------

The geometries and data of the differend adsorbed species are in the folder ``structures``.

    ``Yang_structures``

    This folder contains the xyz files of all adsorbates and transition states included in the reaction model for Rh(111) (``Rh111``) and Rh(211) (``Rh211``) originally published by |Yang2016|_. Only those files are contained, which are a single adsorbate on the respective facet and the site names for the Rh(111) are changed *s* to *t*. The original file names are in the respective xyz comment line.
\

    ``additional_structures``

    This folder contains the xyz files of the adsorbate and transition states additionally calculated on Rh(111) (``Rh111``) and Rh(211) (``Rh211``).

.. |Yang2016| replace:: Yang *et al*
.. _Yang2016: https://doi.org/10.1021/jacs.5b12087

.. |kmos_paper| replace:: kmos
.. _kmos_paper: https://doi.org/10.1016/j.cpc.2014.04.003

.. |kmos| replace:: kmos repository
.. _kmos: https://github.com/m-deimel/kmos.git

.. |kmos_branch| replace:: feature-rate-distribution
.. _kmos_branch: https://github.com/m-deimel/kmos/tree/feature-rate-distribution

.. |CO_hydrogenation| replace:: \"Selectivity Trends and Role of Adsorbate-Adsorbate Interactions in CO Hydrogenation on Rhodium Catalysts\"
.. _CO_hydrogenation: https://doi.org/
