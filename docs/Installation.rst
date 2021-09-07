Installation instructions
***************************

biops can be installed with pip or conda.

.. code-block:: console

   pip install biops

or

.. code-block:: console

   conda install -c omvatten biops

Recommended installation procedure
##################################

The recommended method is to download `Anaconda <https://www.anaconda.com/products/individual>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_

Open the Anaconda or Miniconda prompt (a terminal window).

Create a new environment. You can, for example, call it biops_env. Then, activate the environment.

.. code-block:: console

   conda create -n biops_env python=3.9
   conda activate biops_env

Install biops using:

.. code-block:: console

   conda install pandas numpy matplotlib scipy
   pip install biops

To start using biops, you need some way of writing and executing Python code. I use `Spyder <https://www.spyder-ide.org/>`_ or Jupyter notebooks. You can install Spyder like this:

.. code-block:: console

   conda install spyder

To run Spyder, simply type:

.. code-block:: console

   spyder

You can install Jupyter like this:

.. code-block:: console

   conda install jupyter

To start a Jupyter notebook, simply type:

.. code-block:: console

   jupyter notebook

To check if biops works, you can run the following code:

.. code-block:: python

   import biops
   r = biops.ifas.Reactor()
   r.view_process()
   
Hopefully, this will plot a figure showing a reactor followed by a settler.
