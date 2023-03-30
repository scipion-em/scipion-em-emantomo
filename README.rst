===================================
Scipion plugin for EMAN2 Tomography
===================================

.. image:: https://img.shields.io/pypi/v/scipion-em-emantomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-emantomo
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-emantomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-emantomo
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-emantomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-emantomo
        :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/dm/scipion-em-emantomo
        :target: https://pypi.python.org/pypi/scipion-em-emantomo
        :alt: Downloads

This plugin provide wrappers around several programs of `EMAN2 <https://blake.bcm.edu/emanwiki/EMAN2>`_ tomography software suite.

============
Installation
============
The plugin can be installed in user (stable) or developer (latest, may be unstable) mode:

**1. User (stable) version:**:

.. code-block::

    scipion3 installp -p scipion-em-emantomo

**2. Developer (latest, may be unstable) version:**:

* Clone the source code repository:

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-emantomo.git

* Move to devel branch:

.. code-block::

    git checkout devel

* Install:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-emantomo --devel

=========
Protocols
=========
The integrated protocols are:

1. PCA-K Means classification of subtomograms.

2. Average subtomograms.

3. Subtomograms manual picking.

4. Subtomograms extraction from tomogram.

5. Initial model.

6. Subtomograms refinement.

=====
Tests
=====

The installation can be checked out running some tests. To list all of them, execute:

.. code-block::

     scipion3 tests --grep emantomo

To run all of them, execute:

.. code-block::

     scipion3 tests --grep emantomo --run

To run a specific test, for example, the tests to check the protocol for averaging subtomograms (the following command
can be copied from the test list displayed when listing the tests, as explained above):

.. code-block::

    scipion3 tests emantomo.tests.test_protocols_sta_classic_workflow.TestEmanTomoAverageSubtomogramsStaClassic

========
Tutorial
========
The workflow test generates a emantomo workflow that offers an overview of how to use emantomo.

.. code-block::

    scipion3 tests emantomo.tests.test_protocols_sta_classic_workflow

==========
References
==========

* `Single particle tomography in EMAN2. <https://doi.org/10.1016/j.jsb.2015.04.016>`_
  Jesús G. Galaz-Montoya et al., Journal of Structural Biology, 2015.

* `High resolution single particle refinement in EMAN2.1. <https://doi.org/10.1016/j.ymeth.2016.02.018>`_
  James M. Bell et al., Methods, 2016.


===================
Contact information
===================

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net or open an issue_.

We'll be pleased to help.

*Scipion Team*

.. _issue: https://github.com/scipion-em/scipion-em-emantomo/issues
