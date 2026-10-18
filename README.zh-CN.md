# 球棱柱线积分正演（TLI）

用于球棱柱重力场、磁场正演的 C++17 / MATLAB 实现。普通 TLI 输出势及一至三阶笛卡尔导数共 20 个重力分量；磁场接口输出 10 个量。TLI-MMSI 用于规则经纬网格质量壳与磁化壳，采用近场线积分和远场二维积分。

本包为 **0.1.0-rc.1 发布候选包**，根据 2026-10-18 本地研究目录整理。相关论文作者为 Pan Zhang、Jinsong Du，题目见 [英文首页](README.md)。尚未代为创建远程仓库或发布 Release。

## 快速运行

安装 C++17 编译器和 CMake 3.18 以上版本，在仓库根目录执行：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel 2
ctest --test-dir build -C Release --output-on-failure
```

Visual Studio 下运行 `build/Release/tli_example.exe`；Linux/macOS 通常运行 `build/tli_example`。默认构建和快速测试无需 MATLAB、MPI、OpenMP 或论文数据。

MATLAB 中先添加 `Functions/`，运行 `build_Tesseroid_LineIntegration_mex`，再调用四个包装函数。`P` 后缀使用 MEX 内的 OpenMP，不使用 `parfor`。详见 [构建说明](docs/BUILD.md) 及 [MATLAB 算例](examples/demo_matlab.m)。

## 必须区分的约定

| 接口 | 角度 | 常数与单位 |
| --- | --- | --- |
| `Prism` / `Point` / `estimate_all` | 弧度 | 默认 SI 重力常数；采用米和 kg/m³ 可获得 SI 输出 |
| MATLAB、Work、TLI-MMSI 兼容接口 | 度 | 内部 G=1；旧算例多采用 km，输出需要物理尺度换算 |

局部 x、y、z 分别指北、东、径向向外，势采用正的 Newton 势。两类接口的 20 分量排列也有差别，不能直接按列混用。[API.md](docs/API.md) 列出完整映射和单位换算。

普通 TLI 的均匀球棱柱公式与 TLI-MMSI 的插值、近场近似具有不同适用范围。精确极点另有数值分支；上径向面使用外侧极限，下径向面及侧边界通常拒绝未定义的边界取值。详见 [方法说明](docs/METHOD.md)。

## 论文算例与数据

源码仓库保留论文实验脚本和小型汇总表。另一个数据压缩包应直接解压到仓库根目录，以恢复 `Works_Zhangetal2026/...` 的原有数据路径。

- Test1：高精度 3D/2D/1D 降维对比，需要 Symbolic Math Toolbox，耗时较长。
- Test2：完整球壳验证，C++/MPI 驱动可从模型定义重新生成结果。
- Test3：极轴验证；重新计算解析参考值需要外部 `Tesseroid_AnalyticalPolar_vpa.m`。
- Test4：Laplace/Poisson 验证，数据包含归档参数和结果。
- RunTimeTest1–3：单体及全球重力/磁场耗时实验；全球模型和地图绘图资源的外部依赖已明确列出。

归档结果来自原研究目录，不能自动视作本次发布快照已重新完成全量计算的证据。[复现说明](docs/REPRODUCIBILITY.md)、[验证说明](docs/VALIDATION.md) 和 [发布审阅记录](docs/RELEASE_AUDIT.zh-CN.md) 写明已验证范围与剩余限制。

按用户选择采用 [MIT 许可证](LICENSE)。原目录保持不变；所有打包改动均在发布副本中完成。
