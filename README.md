# 一些小脚本或者程序

* `cgenff_charmm2gmx_py3_nx3.py` 此脚本只支持networkx >= 2.5系列，包括3.x的也可以运行，测试环境Python 3.8和3.10

* `PymixMem.py` 膜构建脚本
【分子动力学磷脂膜构建工具--PymixMem-哔哩哔哩】 https://b23.tv/VuHZFdq

* `nanotube1.5` 纳米管和氧化石墨烯建模工具

* `vmd1.9.3ExtendPlugin_win_20240422.zip`  拓展的Windows版本VMD 1.9.3的分子文件读取库

* 大体系原子近邻查找算法对比
  
  * [MDAnalysis](https://github.com/MDAnalysis/mdanalysis)库本身有三种算法，暴力遍历这里不提，因为大体系根本算不动，周期性kdtree表现比网格划分更加耗时，并且都无法并行计算。
  * [freud](https://freud.readthedocs.io/en/stable/)库本身有两种算法，AABBQuery算法要优于LinkCell算法，分别和mda的网格划分和周期性kdtree速度持平，有TBB并行
  * [Ovito](https://www.ovito.org/)算法主要也是基于网格划分的，性能很好，有并行
  * **cpp**是我整合的一个网格划分算法，基本上也是和Ovito持平，有并行
  * [VMD](http://www.ks.uiuc.edu/Research/vmd/)可视化软件实际上也是基于网格划分算法，性能也不错，不过未考虑周期性边界，此处未列出

小体系并行对计算速度提升不大，甚至并行开销影响速度。


