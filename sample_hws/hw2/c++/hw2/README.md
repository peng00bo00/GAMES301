# hw2

## 要求

**C++17**
需确保 C++ 编译器支持 C++17 标准

**CMake 3.22**

需确保 CMake 版本至少是 3.22

**Gnuplot 5.2.6**
本项目在**运行时**需要使用 Gnuplot 来绘制能量随迭代次数变化图。

如果需要显示该图，请确保 Gnuplot 版本至少是5.2.6。

注意：在编译时并不需要安装该程序。

=== "Ubuntu + GCC"

    ```bash
    sudo apt update
    sudo apt install gnuplot
    ```
    
    !!! note ""
        Or download the latest version from [www.gnuplot.info](http://www.gnuplot.info). If you're using an installer, make sure you mark the option "Add application directory to your PATH environment variable".


=== "Mac Os + Clang"

    ```bash
    brew install gnuplot
    ```
    
    !!! note ""
        Or download the latest version from [www.gnuplot.info](http://www.gnuplot.info). If you're using an installer, make sure you mark the option "Add application directory to your PATH environment variable".

=== "Windows + MSVC"

    !!! warning ""
        Download Gnuplot from [www.gnuplot.info](http://www.gnuplot.info) and install it.
    
        If you're using the Gnuplot installer, make sure you mark the option "Add application directory to your PATH environment variable"

**Library**

本项目编译依赖于以下库：

- libigl
- Eigen
- Matplot++

本项目默认通过 CMake 的 FetchContent 来下载这些依赖库，因此并不需要额外配置。如果想使用本地的库，请确保 Eigen 版本至少为 3.4.0。

## 编译及运行

```bash
mkdir build
cd build
cmake ..
make
./hw2
```

