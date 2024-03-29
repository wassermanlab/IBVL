# IBVL - Open CGA installation

## Software Installation
### Install Python3
#### MacOS
If you don't have Python 3 installed, you can install it on Mac using Homebrew:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
brew install python
```

Check that the installation worked
```
python3 --version
pip3 --version
```

Optionally, you can create an alias to use `python` instead of `python3` by adding the following to your `bashrc`
```
alias python=/usr/local/bin/python3
```
### Clone the repo
If you haven't done so already, clone the repo onto your machine:
```
git clone git@github.com:scorreard/IBVL.git
```

### Set up an Environment
We are using an environment to manage Python libraries to prevent version issues. Two options are to use `virtualenv` or `conda`
#### Using `virtualenv`
```
pip3 install virtualenv
python3 -m venv venv
source venv/bin/activate
```

#### Using Conda
Alternatively, you can install conda [here](https://docs.conda.io/en/latest/miniconda.html) and set up an environment with it
```
conda create -n opencga python=3.8
conda activate opencga
```

## Install requirements
You will need to install the dependencies in `requirements.txt`
```
pip3 install -r requirements.txt
```

## Pull the repository
Pull the repository to have the python scripts
```
git pull
```

## Create the `.env` file
Ask for the `.env` file

A `.env` may look like this :
```
export OPENCGA_HOST="http://192.168.21.116:8080/opencga"
export OPENCGA_USERNAME="username"
export OPENCGA_PASSWORD="password"
```
Source the file
```
source .env
```


---
