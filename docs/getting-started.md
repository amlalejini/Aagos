# Getting started

## One-time setup

You should only need to do these steps once:

1) Clone the github repository on your machine.
   ```
   git clone https://github.com/amlalejini/Aagos.git
   ```
2) Initialize and update the repository's git submodules. `cd` into the Aagos repository that you just cloned, and run
   ```
   git submodule update --init --recursive
   ```

3) Create a Python virtual environment. From inside the repository directory, run
    ```
    python3 -m venv pyenv
    ```
    This will create a [python virtual environment](https://docs.python.org/3/library/venv.html) called `pyenv`.
4) Activate the python virtual environment, and install python dependencies.
    ```
    source pyenv/bin/activate
    ```
    With the virtual environment 'activated', you can install packages/run python programs in the virtual environment.

    The python dependencies for this project are documented in the `requirements.txt` file.
    To install python dependencies, you can run:

    ```
    pip install -r requirements.txt
    ```

## Compiling Aagos

To compile in debug mode:

```
make debug
```

To compile for experimentation:

```
make native
```