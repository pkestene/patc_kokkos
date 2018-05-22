# build the code

We strongly recommend the out-of-source build, so that one can have one build directory per architecture.

```shell
mkdir build
cd build
cmake ..
make
```

You should now have executable *euler2d*. You can run a simply implode test like this
```shell
cd src
./euler2d ./test_implode.ini
```

# Build the code (the old way)

Just for testing, the old Makefile is still there, if by any chance it can be usefull

```shell
make -f Makefile.standalone
```

# Visualize outputs

Output file format is VTK image data (for structured grid).
You can use either *paraview* or directly in python using the script misc/plot_data.py





