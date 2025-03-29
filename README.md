# DNAOAD

To run the examples, enter the `examples` folder and use a terminal. For instance:

```bash
sh create_exec.sh ex1.f90
```

On Windows, use:

```cmd
create_exec.bat ex1.f90
```

Alternatively (and this is the recommended method), you can manually compile and execute the source files.

The examples in the `SourcesF` folder can be compiled directly using the Intel Fortran compiler. For example:

```bash
ifx ex3.f90
```

If you are using **gfortran**, use the `SourcesF_gfortran` folder instead.

For more details, see **Section 4** of the article:  
https://doi.org/10.48550/arXiv.2501.04159
