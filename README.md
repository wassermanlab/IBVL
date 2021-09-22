# IBVL 

## Install Python3
### Install on MacOS
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
brew install python
```

## Setting up an Environment
### Using `virtualenv`
```
pip3 install virtualenv
python3 -m venv venv
source venv/bin/activate
```

### Using Conda
Alternatively, you can install conda and set up an environment with it
```
conda create -n opencga python=3.8
conda activate opencga
```

## Install requirements
```
pip3 install -r requirements.txt
```
