========================
EMAN2 Tomography plugin
========================

This plugin provide wrappers around several programs of `EMAN2 <https://blake.bcm.edu/emanwiki/EMAN2>`_ tomography software suite.

+------------------+------------------+
| stable: |stable| | devel: | |devel| |
+------------------+------------------+

.. |stable| image:: http://scipion-test.cnb.csic.es:9980/badges/eman2_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/eman2_sdevel.svg


Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-emantomo

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-emantomo.git

    * install

    .. code-block::

        scipion installp -p path_to_scipion-em-emantomo --devel

EMAN2 binaries will be installed automatically with the plugin using a Conda environment.

* Default installation path assumed is ``software/em/eman-2.91``, if you want to change it, set *EMANTOMO_HOME* in ``scipion.conf`` file pointing to the folder where the EMANTOMO is installed.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

   scipion3 tests emantomo.tests.test_protocols_emantomo.TestEmanTomoTempMatch
   scipion3 tests emantomo.tests.test_protocols_emantomo.TestEmanTomoSubtomogramRefinement
   scipion3 tests emantomo.tests.test_protocols_emantomo.TestEmanTomoReconstruction
   scipion3 tests emantomo.tests.test_protocols_emantomo.TestEmanTomoInitialModel
   scipion3 tests emantomo.tests.test_protocols_emantomo.TestEmanTomoExtraction

A complete list of tests can also be seen by executing ``scipion test --show --grep emantomo``

Supported versions
------------------

* 2.9
* 2.91 (**Default version**)

Protocols
---------

* `boxer (new interactive e2spt_boxer.py) <https://blake.bcm.edu/emanwiki/EMAN2/Programs/e2tomoboxer>`_
* Template Matching
* Initial model SGD
* Tomogram reconstruction
* Subtomogram extraction
* Subtomogram refinement

References
----------

1. \G. Tang, L. Peng, P.R. Baldwin, D.S. Mann, W. Jiang, I. Rees & S.J. Ludtke. (2007) EMAN2: an extensible image processing suite for electron microscopy. J Struct Biol. 157, 38-46. PMID: 16859925
