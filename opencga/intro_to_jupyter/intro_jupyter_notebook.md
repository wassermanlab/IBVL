# Introduction to Jupyter Notebooks
## What is a Jupyter Notebook?
Jupyter Notebook is a web-based application that integrates code, visualizations, and text into a single document. It allows you to run code, add comments and explanations, view charts and formulas, and display the output of the code all in one place. Although many different programming languages can be used with Jupyter Notebooks, we will demonstrate how to use them with Python. 

More information about Project Jupyter and what Jupyter notebooks are can be found [here](https://jupyter.org/)

More information on how to install and use Jupyter Notebook can be found [here](https://jupyter-notebook.readthedocs.io/en/stable/)


## Installing Jupyter Notebook 
For use with Python, Jupyter Notebook can be installed independently or as part of Anaconda, which is a Python and R distrubution. For this document we will install Jupyter independently using `pip`


### Install Homebrew (MacOS only)
The first step is to install Python3 on your laptop. You can do so using Homebrew, which is a package management system for MacOS. In order to install brew, you will need to open a new terminal window. 

Navigate to either the spotlight search or launchpad and type in "terminal". Click on the terminal icon to open a new window. The terminal is program that allows you to interact and perform tasks on the computer without using a graphical interface. We will use the terminal to run our Jupyter notebooks. 

After opening a new terminal window, run the following command to install Homebrew:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```


### Install Python 3 (MacOS only)
Once you have installed Homebrew, you can now install Python and pip:
```
brew install python
```


### Install Python 3 (Ubuntu only)
Open a new terminal window and run the following command:
```
sudo apt update
sudo apt-get install python3-pip
```


### Verify Python Installation Worked
You can check that the installation worked by running the following command:
```
python3 --version
pip3 --version
```

The outputs of the above command should be:
```
Python 3.9.5
pip 21.1.3
```

It is okay if you have different version numbers, but they should be equal to or greater than the above numbers. 


### Install Jupyter
Now that you have Python installed, navigate to the folder containing your Jupyter notebook. You can change directories by using the `cd` command followed by the path to your folder:
```
cd ~/Documents/IBVL/opencga
```

You will need to create a virtual environment for Python. This ensures that the libraries you install for one project won't interfere with the libraries you install for another project by introducing version issues. 

Install `virtualenv` to manage your virtual environment:
```
pip3 install virtualenv
```

Once you have installed `virtualenv`, create a new virtual environment:
```
python3 -m venv venv
```

The above command created a new virtual environment in your current folder called `venv`. You can activate this environment by navigating to the same directory for which it is located and running the following:
```
source venv/bin/activate
```

Now that you have activate your virtual environment, you can install Jupyter notebook
```
pip3 install notebook
```

You can run Jupyter Notebook with the following command:
```
jupyter notebook
```

This should open a window in your internet browser where you are able to interact with your notebook!


## Using Jupyter Notebooks
When you want to start using a Jupyter notebook, open a terminal window run the following command:
```
jupyter notebook
```

Once the web browser opens with Jupyter, you can click on folders until you find the notebook you want to work with. Jupyter Notebooks will have the extension `.ipynb`. Clicking on one of these files will open the notebook in your web browser. For this part of the tutorial we will refer to the `example_notebook.ipynb` notebook. 

Notebooks are made up of sequential cells, with each cell being a multiline text input field. You can run a cell using `Ctrl-Enter` on your keyboard, or you can click the "Play" button at the top of the notebook. 

There are three types of cells in a notebook:
1. code cells
2. markdown cells
3. raw cells

The default cell type is a code cell, but you can change the cell type using the toolbar under `Cell>Cell Type`.

You can insert a new cell below your current cell using `Shift-Enter`.

### Code Cells
Code cells allow you to edit and write blocks of code, similar to your text editor. When you run a code cell the code in that cell is executed and the results are displayed in the notebook. 

The cell output can be text, figures, tables, charts, and much more. 

The first code cell in our example notebook is a line of Python code which prints "hello world!" to the output. If you run this cell using `Ctrl-Enter`, you will see "hello world!" displayed underneath the code cell.  

### Markdown Cells
Markdown cells provide a way for you to include text in your document using the Markdown language. This allows you to perform text markup on the text within your cells, such as bolding text, creating lists, adding links, etc. For more information on the syntax involved in markdown, refer to [this cheat sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet). 

The first cell in the example notebook is a markdown cell. When you double click on a markdown cell you will notice the text changes. This allows you to edit the text within the cell. Once you are done editing the text, run `Ctrl-Enter` to execute the changes.

### Raw Cells
Raw cells are not used very often, as they are not evaluated by the notebook. You can use raw cells if you intend to use the `nbconvert` tool to convert your notebook to another format (such as HTML or Latex). When this happens, raw cells will be converted in a way specific to the output you want. 


## Closing Jupyter Notebook
When you are done using Jupyter Notebook, save your changes and exit out of the web browser. You can then navigate to your terminal and run `Ctrl-c`. You will be asked "Shutdown this notebook server (y/[n])?", to which you can press `y-Enter`.

You can now exit out of the terminal window.