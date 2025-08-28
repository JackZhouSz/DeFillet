<p align="center">
  <p align="center">
    <h1 align="center">DeFillet: Detection and Removal of Fillet Regions in Polygonal CAD Models</h1>
  </p>


<p align="center">
  <a href="https://raw.githubusercontent.com/xiaowuga/xiaowuga.github.io/main/pub/static/pdf/Defillet_Sig_2025.pdf">
    <img src="https://img.shields.io/badge/Paper-PDF-red?logo=adobeacrobatreader&logoColor=white" alt="Paper">
  </a>
  <a href="https://xiaowuga.github.io/pub/DeFillet.html">
    <img src="https://img.shields.io/badge/Project_Page-Website-green?logo=googlechrome&logoColor=white" alt="Project Page">
  </a>
  <a href="https://www.bilibili.com/video/BV1aCtMzPEXe/?spm_id_from=333.1007.top_right_bar_window_history.content.click&vd_source=092295aa747638ab207808257f039dea">
    <img src="https://img.shields.io/badge/%F0%9F%8E%A5%20Video-Bilibili-blue?logo=bilibili&logoColor=white" alt="Video">
  </a>

[//]: # (  <a href="https://github.com/xiaowuga/DeFillet">)

[//]: # (<img src="https://img.shields.io/badge/Poster-PDF-purple?logo=adobeacrobatreader&logoColor=white" alt="Poster">)

[//]: # (  </a>)
</p>



<p align="center" style="font-size:16px">
    <a target="_blank" href="https://xiaowuga.github.io/"><strong>Jing-En Jiang</strong></a>
    ·
    <a target="_blank" href="https://github.com/hanxiaowang00"><strong>Hanxiao Wang</strong></a>
    ·
    <a target="_blank" href="https://zikai1.github.io/"><strong>Mingyang Zhao</strong></a>
    ·
    <a target="_blank" href="https://sites.google.com/site/yandongming/"><strong>Dong-Ming Yan</strong></a>
    ·
    <a target="_blank" href="https://xk.qust.edu.cn/info/1041/4695.htm"><strong>Shuangmin Chen</strong></a>
    .
    <a target="_blank" href="https://irc.cs.sdu.edu.cn/~shiqing/index.html"><strong>Shiqing Xin</strong></a>
    .
    <a target="_blank" href="https://faculty.sdu.edu.cn/tuzhanghe/en/index.htm"><strong>Changhe Tu</strong></a>
    .
    <a target="_blank" href="https://engineering.tamu.edu/cse/profiles/Wang-Wenping.html"><strong>Wenping Wang</strong></a>
</p>

![](./asset/teaser.gif)








This repository contains the official implementation of our SIGGRAPH 2025 paper "DeFillet: Detection and Removal of Fillet Regions in Polygonal CAD Models".


**Please give a star and cite if you find this repo useful.**

## Platform
- Windows 11
- CLion2024.1.2 +  Visual Stdio 2022
- Intel(R) Core i9-13900K

## Dependence

The dependent libraries of our code includes:
- Eigen3 (3.4.0 or later): Solving linear equations.
- libigl (2.5.0 or later): Solving linear equations.
- Easy3D (2.5.2 or later): IO operation and data structure for mesh.
- CGAL (5.6 or later): Voronoi diagram computation.
- CLI11 (2.4.0 or later): Command line.
- nlohmann_json (3.11.3 or later): Read json file.

For Eigen3, libigl, CGAL, CLI11 and nlohmann_json, we recommend using [vcpkg](https://github.com/microsoft/vcpkg) to install them.
```shell
# Eigen3
vcpkg install Eigen3:x64-windwos
# libigl
vcpkg install libigl:x64-windwos
# CGAL
vcpkg install cgal:x64-windwos
# cli11
vcpkg install cli11:x64-windows
# nlohmann_json
vcpkg install nlohmann_json:x64-windows
```

For Easy3D, please visit the [Prof.Nan's repository](https://github.com/LiangliangNan/Easy3D) to obtain it and follow the instructions for installation. You need to set the `${Easy3D_DIR}` environment variable to the directory that contains the `Easy3DConfig.cmake` file.


## How to Build

Building our code in CLion:
```
# File -> Setting -> Build, Execution, Deployment -> CMake -> CMake Option :
-DCMAKE_TOOLCHAIN_FILE=${YOUR_VCPKG_INSTALL_PATH}/scripts/buildsystems/vcpkg.cmake
```
Making sure that your following settings are correct:
- Toolchains : `Visual Stdio`
- Architecture : `amd64`
- Build Type : `release`

## Usage

### 1.detector_cli

```shell
detector_cli.exe -c detector_config.json
```

The executable target `detector_cli.exe` takes a `.ply` mesh file as input and outputs the fillet segmentation results and intermediate results to a specified output folder.
You can modify the `detector_config.json` file to change the input and output folder paths as well as algorithm parameters (configured based on the details in the paper).

### 2.removal_cli

Coming soon!

## Citation
If you make use of our work, please cite our paper:

```bibtex
@article{jiang2025defillet,
    author = {Jiang, Jing-En and Wang, Hanxiao and Zhao, Mingyang and Yan, Dong-Ming and Chen, Shuangmin and Xin, Shiqing and Tu, Changhe and Wang, Wenping},
    title = {DeFillet: Detection and Removal of Fillet Regions in Polygonal CAD Models},
    journal={ACM Transactions on Graphics (TOG)},
    publisher = {Association for Computing Machinery},
    year = {2025},
    address = {New York, NY, USA},
    volume = {44},
    number = {4},
    issn = {0730-0301},
    url = {https://doi.org/10.1145/3731166},
    doi = {10.1145/3731166} 
}
```


## Maintaince

If any problem, please contact me via <xiaowuga@gmail.com>.


## License
DeFillet is under AGPL-3.0, so any downstream solution and products (including cloud services) that include DeFillet code inside it should be open-sourced to comply with the AGPL conditions. For learning purposes only and not for commercial use. If you want to use it for commercial purposes, please contact us first.





