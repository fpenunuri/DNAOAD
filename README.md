# DNAOAD

**DNAOAD** (Dual Number Arbitrary Order Automatic Differentiation) is a Fortran implementation of dual numbers for arbitrary order automatic differentiation.

## 🧪 Running Examples

To run the examples, navigate to the `examples` directory and use a terminal. For example:

```sh
sh create_exec.sh ex1.f90
```

On **Windows**, use:

```bat
create_exec.bat ex1.f90
```

> **Note:** Alternatively (and this is the recommended method), you can manually compile and run the source files.

## ⚙️ Manual Compilation

The examples in the `SourcesF` folder can be compiled directly using the **Intel Fortran** compiler:

```sh
ifx ex3.f90
```

If you're using **gfortran**, use the `SourcesF_gfortran` folder instead.

## 📖 More Information

For further details, see **Section 4** of the article:

📎 [https://doi.org/10.48550/arXiv.2501.04159](https://doi.org/10.48550/arXiv.2501.04159)

## 🔗 Related repository

See also: [DNAOAD-FPM](https://github.com/fpenunuri/DNAOAD-FPM), which provides an implementation compatible with the Fortran Package Manager (FPM).


