Installation instructions
***************************
Download `Anaconda <https://www.anaconda.com/products/individual>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_

Open the Anaconda or Miniconda prompt (a terminal window).

Create a new environment. You can, for example, call it biops_env. Then, activate the environment.

.. code-block:: console

   conda create -n biops_env
   conda activate biops_env

biops depends on four other python packages: pandas, numpy, matplotlib, and scipy. To install them, use the following command:

.. code-block:: console

   conda install pandas numpy matplotlib scipy
   
Next, install biops using pip.

.. code-block:: console

   pip install biops

To start using biops, you need some way of writing and executing Python code. I use `Spyder <https://www.spyder-ide.org/>`_. You can install Spyder like this:

.. code-block:: console

   conda install spyder

To run Spyder, simply type:

.. code-block:: console

   spyder

To check if biops works, you can run the following code:

.. code-block:: python

   import biops
   r = biops.ifas.Reactor()
   r.view_process()
   
Hopefully, this will plot a figure showing a reactor followed by a settler.